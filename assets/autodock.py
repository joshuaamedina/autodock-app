from vina import Vina
from mpi4py import MPI
import subprocess
import pickle
import os
import blosc
from os.path import exists
from os.path import basename
from sys import argv
import argparse

parser = argparse.ArgumentParser(description ='vina commands')
parser.add_argument('-r','--receptor',type=str, help='Receptor for processing')
parser.add_argument('-c','--center',type=str,required=True,help='center position')
parser.add_argument('-s','--size',type=str,required=True,help='box size')
parser.add_argument('-m','--module',type=str,required=True,help='module for docking')
parser.add_argument('-d','--docking',type=str,required=True,help='basic or flexible docking')

args = parser.parse_args()
docking = args.docking
center = (args.center).split(",")
center_x = float(center[0])
center_y = float(center[1])
center_z = float(center[2])
box = (args.size).split(",")
size_x = float(box[0])
size_y = float(box[1])
size_z = float(box[2])

# Setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

receptor=args.receptor
flex_receptor=f'{receptor}_flex'
receptor_path = f'/scratch/09252/adarrow/autodock/input/receptors'
if docking == 'basic':
    flexible = False
else:
    flexible = True
sidechains = ['THR315']
docking_type = args.module
ligand_library = '/scratch/02875/docking/test/Enamine-PC/100_sets'
config_path = './configs/config.config'
user_configs = {'center_x': center_x, 'center_y': center_y, \
                'center_z': center_z, 'size_x': size_x, \
                'size_y': size_y, 'size_z': size_z}
'''
def prep_config():
    with open(config_path, 'w+') as f:
        for config, value in user_configs.items():
            f.write(f'{config} = {value}\n')
    f.close()
'''
def prep_maps():
    if docking_type == 'ad4':
        if exists(f'{receptor}.gpf'):
            subprocess.run([f"rm {receptor}.gpf"], shell=True)
        subprocess.run([f"python3 ./scripts/write-gpf.py --box {config_path} {receptor_path}/{receptor}.pdbqt"], shell=True)
        subprocess.run([f"./scripts/autogrid4 -p {receptor}.gpf"], shell=True)

def prep_receptor():
    if exists(f'{receptor_path}/{receptor}H.pdb'):
        subprocess.run([f'prepare_receptor -r {receptor_path}/{receptor}H.pdb -o {receptor_path}/{receptor}.pdbqt'], shell=True)
    if flexible == True:
        subprocess.run([f"pythonsh prepare_flexreceptor.py -g {receptor}.pdbqt -r {receptor_path}/{receptor}.pdbqt -s {'_'.join(sidechains)}"], shell=True)
        subprocess.run([f"mv *receptor* {receptor_path}"], shell=True)

def prep_ligands():
    # Returns a list where each item is the path to a pickled and compressed text file containing multiple ligand strings
    ligand_paths = []
    for dirpath, dirnames, filenames in os.walk(ligand_library):
        for filename in filenames:
            ligand_paths.append(f'{dirpath}/{filename}')
    return ligand_paths

def run_docking(ligands, v):
    for index, filename in enumerate(ligands):
        ligand = ligands[filename]
        v.set_ligand_from_string(ligand)
        v.dock()
        v.write_poses(f'output_{filename}', n_poses=1, overwrite=True)
        subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' output_{filename} \
                        | awk '{{print $4}}' >> results_{rank}.txt; echo {filename} \
                        >> results_{rank}.txt"], shell=True)

def unpickle_and_decompress(path_to_file):
    with open(path_to_file, 'rb') as f:
        compressed_pickle = f.read()
    depressed_pickle = blosc.decompress(compressed_pickle)
    dictionary_of_ligands = pickle.loads(depressed_pickle)
    return dictionary_of_ligands

def pre_processing():
    #prep_config()
    prep_receptor()
    prep_maps()
    ligands = prep_ligands()
    return ligands

def processing():
    # Initialize Vina or AD4 configurations
    if docking_type == 'vina':
        v = Vina(sf_name='vina', cpu=4, verbosity=0)
        if flexible == True:
            v.set_receptor(f'{receptor}.pdbqt', f'{flex_receptor}.pdbqt')
        else:
            v.set_receptor(f'{receptor}.pdbqt')
        uc = user_configs
        v.compute_vina_maps(center=[float(uc['center_x']), float(uc['center_y']), float(uc['center_z'])], \
                            box_size=[float(uc['size_x']), float(uc['size_y']), float(uc['size_z'])])
    elif docking_type == 'ad4':
        v = Vina(sf_name='ad4', cpu=0, verbosity=0)
        v.load_maps(map_prefix_filename = receptor)
        
    # Ask rank 0 for ligands and dock until rank 0 says done
    while True:
        comm.send(rank,dest=0) # Ask rank 0 for another set of ligands
        ligand_set_path = comm.recv(source=0) # Wait for a response
        if ligand_set_path == 'no more ligands':
            comm.send('message received--proceed to post-processing',dest=0)
            break
        # Pickle load ligand set
        ligands = unpickle_and_decompress(ligand_set_path)
        # Dock each ligand in the set
        run_docking(ligands, v)



def sort():
    subprocess.run(["cat results* >> results_merged.txt"], shell=True)
    INPUTFILE = 'results_merged.txt'
    OUTPUTFILE = 'processed_results.txt'
    
    result = []

    with open(INPUTFILE) as data:
        line = data.readline()
        while line:
            filename = basename(line.split()[-1])
            v = data.readline().split()[0]
            result.append(f'{v} {filename}\n')
            line = data.readline()

    with open(OUTPUTFILE, 'w') as data:
        data.writelines(sorted(result, key=lambda x: float(x.split()[1])))
    
    subprocess.run(["rm results*; mv *map* *.gpf ./output/maps"], shell=True)

def main():
    if rank == 0:
        # Pre-Processing
        ligands = pre_processing()
        # Let other ranks know pre-processing is finished; they can now ask for work
        for i in range(size):
            comm.sendrecv('pre-processing finished; ask for work', dest=i)

        # Until all ligands have been docked, send more work to worker ranks
        count = 0
        while count < 1:
            source = comm.recv(source=MPI.ANY_SOURCE)
            comm.send(ligands.pop(), dest=source)
            count+=1
#        while ligands:
#            source = comm.recv(source=MPI.ANY_SOURCE)
#            comm.send(ligands.pop(), dest=source)

        # When all ligands have been sent, let worker ranks know they can stop
        for i in range(size):
            comm.send('no more ligands', dest=i)
            comm.recv(source=i)

        # Post-Processing
        sort()
    else: # All ranks besides rank 0
        comm.recv(source=0) # Wait for rank 0 to finish pre-processing
        comm.send(rank, dest=0)
        processing()    
        


main()
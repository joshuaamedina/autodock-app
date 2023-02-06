from vina import Vina
from mpi4py import MPI
import subprocess
import pickle # For unpickling ligand files
import os
import blosc # For decompressing compressed ligand files
from os.path import exists
from os.path import basename # Used in the sorting function
import argparse # To accept user inputs as command line arguments

"""
Setup base MPI declarations
In this MPI implementation, rank 0 acts as the director that all other ranks 
work for. Rank 0 first completes all pre-processing. As a final step, it 
creates a list containing full filepaths to all ligand files in the ligand 
library. Then it accepts messages from worker ranks saying they are ready for 
more work. It sends work for as long as there are ligand files left. Then rank 
0 waits to hear from all ranks that they have finished processing, then 
proceeds to do the post-processing work. 
"""
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Assign user inputs
parser = argparse.ArgumentParser(description ='vina commands')
parser.add_argument('-r','--receptor',type=str,
                    help='Receptor for processing')
parser.add_argument('-c','--center',type=str,required=True,
                    help='center position')
parser.add_argument('-s','--size',type=str,required=True,
                    help='box size')
parser.add_argument('-m','--module',type=str,required=True,
                    help='module for docking')
parser.add_argument('-d','--docking',type=str,required=True,
                    help='basic or flexible docking')
parser.add_argument('-f','--sidechains',type=str,required=False,
                    help='sidechains if flexible docking')
parser.add_argument('-n','--number',type=int,required=True,
                    help='top n results to return')
parser.add_argument('-ll','--ligand_library',type=str,required=True,
                    help='ligand library')

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

full_receptor=args.receptor
receptor=full_receptor.split('.')[0]
flex_receptor=f'{receptor}_flex'
if docking == 'basic':
    flexible = False
elif docking == 'flexible':
    flexible = True
    sidechains = (args.sidechains).split()
docking_type = args.module
ligand_library = args.ligand_library
config_path = './configs/config.config'
user_configs = {'center_x': center_x,
                'center_y': center_y,
                'center_z': center_z,
                'size_x': size_x,
                'size_y': size_y,
                'size_z': size_z}

number_of_outputs = args.number

# Internal variables
cpus = 4
verbosity = 0
poses = 1

def check_user_configs():
    assert (size_x <= 30 and size_y <= 30 and size_z <=30), \
           "box size is outside the bounds (30)" 
    assert (full_receptor.endswith('.pdb') or full_receptor.endswith('.pdbqt')), \
           "Please provide a .pdb or .pdbqt file"

def check_user_sidechains():
    # Check that user-inputted flexible sidechains are found within the .pdbqt file
    if flexible == False:
        return
    all_sidechains = []
    with open(f'{receptor}.pdbqt', 'r') as r:
        line = r.readline()
        while line:
            line = r.readline()
            if line.startswith('ATOM'):
                 all_sidechains.append(line.split()[3] + line.split()[5])
    for sidechain in sidechains:
        assert sidechain in all_sidechains, \
               "Please provide valid flexible sidechain names (e.g. THR315)"


def prep_config():
    # Write user inputted configurations (grid center and box size) to a 
    #   config file for AutoDock Vina to use
    with open(config_path, 'w+') as f:
        for config, value in user_configs.items():
            f.write(f'{config} = {value}\n')
    f.close()


def prep_maps():
    # Generate affinity maps using write-gpf.py and AFDRSuite's autogrid4. 
    #   Generates maps for all possible ligand atom types for a given receptor
    if docking_type == 'ad4':
        if exists(f'{receptor}.gpf'):
            subprocess.run([f"rm {receptor}.gpf"], shell=True)
        subprocess.run([f"python3 ./scripts/write-gpf.py --box {config_path} \
                        {receptor}.pdbqt"], shell=True)
        subprocess.run([f"autogrid4 -p {receptor}.gpf"], shell=True)


def prep_receptor():
    # Converts a PDB receptor to a PDBQT, if needed. If the user has specified 
    #   flexible docking, also prepares the rigid receptor and user-chosen 
    #   flexible sidechains.
    if exists(f'{receptor}.pdb'):
        subprocess.run([f'prepare_receptor -r {receptor}.pdb \
                        -o {receptor}.pdbqt'], shell=True)
    if flexible == True:
        subprocess.run([f"pythonsh ./scripts/prepare_flexreceptor.py \
                        -g {receptor}.pdbqt -r {receptor}.pdbqt \
                        -s {'_'.join(sidechains)}"], shell=True)


def prep_ligands():
    # Returns a list where each item is the path to a pickled and compressed 
    #   text file containing multiple ligand strings
    ligand_paths = []
    for dirpath, dirnames, filenames in os.walk(ligand_library):
        for filename in filenames:
            ligand_paths.append(f'{dirpath}/{filename}')
    return ligand_paths


def run_docking(ligands, v, directory):
    # Runs AutoDock on each ligand in the given set; outputs a .pdbqt file 
    #   showing the pose and all scores; appends the ligand name (filename) 
    #   and it's best pose/score to a temporary results file
    output_directory = f'./output/pdbqt/{rank}{directory}'
    if not exists(output_directory):
        os.makedirs(output_directory)
    for index, filename in enumerate(ligands):
        ligand = ligands[filename]
        v.set_ligand_from_string(ligand)
        v.dock()
        v.write_poses(f'{output_directory}/output_{filename}', \
                      n_poses=poses, overwrite=True)
        subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' \
                        {output_directory}/output_{filename} \
                        | awk '{{print $4}}' >> results_{rank}.txt; echo {filename} \
                        >> results_{rank}.txt"], shell=True)


def unpickle_and_decompress(path_to_file):
    # Given a filepath, decompresses and unpickles the file. Returns the 
    #   contents of the file as a dictionary where keys are ligand filenames 
    #   and values are the actual ligand strings
    with open(path_to_file, 'rb') as f:
        compressed_pickle = f.read()
    depressed_pickle = blosc.decompress(compressed_pickle)
    dictionary_of_ligands = pickle.loads(depressed_pickle)
    return dictionary_of_ligands


def pre_processing():
    # Helper function to reduce clutter in main()
    subprocess.run(['mkdir -p configs output/pdbqt output/top_results'], shell=True)
    check_user_configs()
    prep_config()
    prep_receptor()
    prep_maps()
    check_user_sidechains()


def processing():
    # Long method; difficult to de-couple any one part without breaking many things

    # Initialize docking configurations
    # Note: Internal benchmarks showed that on a given node with 128 cores, 
    #   setting ibrun -np to 32 and, below, setting CPU=4 granted the fastest 
    #   docking. Verbosity set to 0 to increase speeds, as stdout cannot be 
    #   captured from Vina's Python module.
    if docking_type == 'vina':
        v = Vina(sf_name='vina', cpu=cpus, verbosity=verbosity)
        if flexible == True:
            v.set_receptor(f'{receptor}.pdbqt', f'{flex_receptor}.pdbqt')
        else:
            v.set_receptor(f'{receptor}.pdbqt')
        uc = user_configs
        v.compute_vina_maps(center=[float(uc['center_x']),
                                    float(uc['center_y']),
                                    float(uc['center_z'])],
                            box_size=[float(uc['size_x']),
                                      float(uc['size_y']),
                                      float(uc['size_z'])])
    elif docking_type == 'ad4':
        v = Vina(sf_name='ad4', cpu=cpus, verbosity=verbosity)
        v.load_maps(map_prefix_filename = receptor)
        
    # Ask rank 0 for ligands and dock until rank 0 says done
    count = 1
    directory = 1
    while True:
        comm.send(rank,dest=0) # Ask rank 0 for another set of ligands
        ligand_set_path = comm.recv(source=0) # Wait for a response
        if ligand_set_path == 'no more ligands':
            comm.send('message received--proceed to post-processing',dest=0)
            break
        ligands = unpickle_and_decompress(ligand_set_path)
        run_docking(ligands, v, directory)
        count += 1
        if count == 100:
            count = 1
            directory += 1



def sort():
    # Cats all results files into one, arranges each line to read: 
    #   (ligand, top score), then sorts by score so that highest scoring 
    #   ligands are on top; prints these sorted results are written to 
    #   processed_results.txt; finally cleans up the directory
    subprocess.run(["cat results* >> results_merged.txt"], shell=True)
    INPUTFILE = 'results_merged.txt'
    OUTPUTFILE = './output/top_results/processed_results.txt'
    
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
    

def isolate_output():
    # Copies the user-specified top n ligand output files to a single directory
    top_ligand_filenames = []
    
    with open('./output/top_results/processed_results.txt', 'r') as results:
        for index, line in enumerate(results):
            if index == (number_of_outputs): 
                break
            top_ligand_filenames.append(line.split()[0])

    for dirpath, dirnames, filenames in os.walk('./output/pdbqt'):
        for top_filename in top_ligand_filenames:
            for filename in filenames:
                if filename == f'output_{top_filename}':
                    subprocess.run([f'mv {dirpath}/{filename} \
                                    ./output/top_results'], shell=True)
    subprocess.run([f'tar -czf top_results.tar.gz ./output/top_results'],
                   shell=True)


def reset():
    # Reset directory by removing all created files from this job
    for dirpath, dirnames, filenames in os.walk('.'):
        for filename in filenames:
            if filename.endswith(('.map', '.txt', '.gpf', '.fld', '.xyz')): 
                subprocess.run([f'rm {dirpath}/{filename}'], shell=True)
    subprocess.run(['rm -r ./output/ ./configs'], shell=True)


def main():
    if rank == 0:
        # Pre-Processing
        pre_processing()
        ligands = prep_ligands()
        # Let other ranks know pre-processing is finished; they can now ask for work
        for i in range(1,size):
            comm.sendrecv('pre-processing finished; ask for work', dest=i)

        # Until all ligands have been docked, send more work to worker ranks
        while ligands:
            source = comm.recv(source=MPI.ANY_SOURCE)
            comm.send(ligands.pop(), dest=source)

        # When all ligands have been sent, let worker ranks know they can stop
        for i in range(1,size):
            comm.send('no more ligands', dest=i)
            comm.recv(source=i)

        # Post-Processing
        sort()
        isolate_output()
        reset()

    else: # All ranks besides rank 0
        comm.recv(source=0) # Wait for rank 0 to finish pre-processing
        comm.send(rank, dest=0)
        processing()    
        


main()

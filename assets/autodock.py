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
if docking == 'rigid':
    flexible = False
elif docking == 'flexible':
    flexible = True
    sidechains = (args.sidechains).split('_')
docking_type = args.module
ligand_library = args.ligand_library
config_path = './configs/config.config'
user_configs = {'center_x': center_x,
                'center_y': center_y,
                'center_z': center_z,
                'size_x': size_x,
                'size_y': size_y,
                'size_z': size_z}


number_of_outputs = args.number if args.number <= 1000 else 1000

# Internal variables
cpus = 4
verbosity = 0
poses = 1
exhaustiveness = 8

def check_user_configs():
    # Check that user-inputted flexible sidechains are found within the .pdbqt file
    for size in [size_x, size_y, size_z]:
        if not (size <= 30 and size >= 1):
           subprocess.run(["echo 'box size is outside the bounds (1-30)' \
                           >> error.txt"], shell=True)
           comm.Abort()
    if not (full_receptor.endswith('.pdb') or full_receptor.endswith('.pdbqt')):
        subprocess.run(["echo 'Please provide a .pdb or .pdbqt file' \
                        >> error.txt"], shell=True)
        comm.Abort()
          
    all_sidechains = []
    xbounds = []
    ybounds = []
    zbounds = []
    with open(f'{receptor}.pdbqt', 'r') as r:
        line = r.readline()
        while line:
            line = r.readline()
            if line.startswith('ATOM'):
                xbounds.append(float(line.split()[6]))
                ybounds.append(float(line.split()[7]))
                zbounds.append(float(line.split()[8]))
                all_sidechains.append(line.split()[3] + line.split()[5])
    if flexible == True:
        for sidechain in sidechains:
            if not sidechain in all_sidechains:
                subprocess.run(["echo 'Please provide valid flexible sidechain \
                                names, separated by spaces (e.g. THR315 GLU268)' \
                                >> error.txt"], shell=True)
                comm.Abort()
                  
    if not (min(xbounds)) <= center_x <= (max(xbounds)):
        subprocess.run(["echo 'Center x coordinate is not within bounds' \
                        >> error.txt"], shell=True)
        comm.Abort()
    if not (min(ybounds)) <= center_y <= (max(ybounds)):
        subprocess.run(["echo 'Center y coordinate is not within bounds' \
                        >> error.txt"], shell=True)
        comm.Abort()
    if not (min(zbounds)) <= center_z <= (max(zbounds)):
        subprocess.run(["echo 'Center z coordinate is not within bounds' \
                        >> error.txt"], shell=True)
        comm.Abort()


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
        try:
            subprocess.run([f'prepare_receptor -r {receptor}.pdb \
                            -o {receptor}.pdbqt'], shell=True)
        except:
            subprocess.run([f"echo 'error on rank {rank}: error prepping receptor' \
                            >> errors.txt"],shell=True)
            comm.Abort()
        
    if flexible == True:
        try:
            subprocess.run([f"pythonsh ./scripts/prepare_flexreceptor.py \
                            -g {receptor}.pdbqt -r {receptor}.pdbqt \
                            -s {'_'.join(sidechains)}"], shell=True)
        except:
            subprocess.run([f"echo 'error on rank {rank}: error prepping receptor' \
                            >> errors.txt"],shell=True)
            comm.Abort()

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
        v.dock(exhaustiveness=exhaustiveness)
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
    subprocess.run(['mkdir -p configs output/pdbqt output/results/ligands'], shell=True)
    prep_config()
    prep_receptor()
    prep_maps()
    check_user_configs()


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

        try:
            ligands = unpickle_and_decompress(ligand_set_path)
        except:
            subprocess.run([f"echo 'error on rank {rank}: could not \
                            unpickle/decompress {ligand_set_path}' \
                            >> errors.txt"],shell=True)
        try:
            run_docking(ligands, v, directory)
        except:
            subprocess.run([f"echo 'error on rank {rank}: docking error with \
                            ligand set {ligand_set_path}, ligands {ligands}' \
                            >> errors.txt"],shell=True)
        count += 1
        if count == 100:
            count = 1
            directory += 1



def sort():
    # Cats all results files into one, arranges each line to read: 
    #   (ligand, top score), then sorts by score so that highest scoring 
    #   ligands are on top; prints these sorted results are written to 
    #   sorted_scores.txt; finally cleans up the directory
    subprocess.run(["cat results* >> results_merged.txt"], shell=True)
    INPUTFILE = 'results_merged.txt'
    OUTPUTFILE = './output/results/sorted_scores.txt'
    
    result = []

    with open(INPUTFILE) as data:
        line = data.readline()
        while line:
            filename = basename(line.split()[-1])
            v = data.readline().split()[0]
            result.append(f'{v} {filename}\n')
            line = data.readline()

    with open(OUTPUTFILE, 'w') as data:
        data.writelines(sorted(result[:number_of_outputs], \
                        key=lambda x: float(x.split()[1])))
    

def isolate_output():
    # Copies the user-specified top n ligand output files to a single directory
    top_ligand_filenames = []
    
    with open('./output/results/sorted_scores.txt', 'r') as results:
        for index, line in enumerate(results):
            top_ligand_filenames.append(line.split()[0])

    for dirpath, dirnames, filenames in os.walk('./output/pdbqt'):
        for top_filename in top_ligand_filenames:
            for filename in filenames:
                if filename == f'output_{top_filename}':
                    subprocess.run([f'mv {dirpath}/{filename} \
                                    ./output/results/ligands'], \
                                    shell=True)
    
    with open('./output/results/combined_docked_ligands', 'w+') as combined:
        for dirpath, dirnames, filenames in os.walk('./output/results/ligands'):
            for filename in filenames:
                lines = open(f'{dirpath}/{filename}', 'r').read()
                combined.write(str(lines))
    
    subprocess.run([f'tar -czf results.tar.gz ./output/results'],
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
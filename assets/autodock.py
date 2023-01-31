from vina import Vina
from mpi4py import MPI
import subprocess
import pickle
import time
import os
import sys
from os.path import exists
from os.path import basename

# Setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
 
receptor='1iep_receptor'
docking_type = 'vina'
ligand_library = ''

def prep_maps(receptor):
    if docking_type == 'ad4':
        if exists(f'{receptor}.gpf'):
            subprocess.run([f"rm {receptor}.gpf"], shell=True)
        subprocess.run([f"python3 ./scripts/write-gpf.py --box ./configs/config.config ./input/receptors/{receptor}.pdbqt"], shell=True)
        subprocess.run([f"./scripts/autogrid4 -p {receptor}.gpf"], shell=True)

def prep_receptor(receptor):
    if exists(f'{receptor}H.pdb'):
        subprocess.run([f'prepare_receptor -r {receptor}H.pdb -o {receptor}.pdbqt'], shell=True)

def run_docking(ligand, v, filename):
    v.set_ligand_from_string(ligand)
    v.dock()
    v.write_poses(f'output_{filename}', n_poses=1, overwrite=True)

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
    # Pre-Processing (only rank 0)
    ligands = None
    if rank == 0:
        prep_receptor(receptor)
        prep_maps(receptor)
        ligands = pickle.load(open('ligands_10.pkl', 'rb'))
    # Initialize Vina or AD4 configurations
    ligands = comm.bcast(ligands, root=0)

    if docking_type == 'vina':
        v = Vina(sf_name='vina', cpu=4, verbosity=0)
        v.set_receptor(f'{receptor}.pdbqt')
        v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
    elif docking_type == 'ad4':
        v = Vina(sf_name='ad4', cpu=0, verbosity=0)
        v.load_maps(map_prefix_filename='1iep_receptor')

    # Run docking
    for index, filename in enumerate(ligands):
        ligand = ligands[filename]
        if(index % size == rank):
            run_docking(ligand, v, filename)
            subprocess.run([f"grep -i -m 1 'REMARK VINA RESULT:' output_{filename} \
                            | awk '{{print $4}}' >> results_{rank}.txt; echo {filename} \
                            >> results_{rank}.txt"], shell=True)
    comm.Barrier()
    # Post-Processing (only rank 0)
    if rank == 0:
        sort()

main()

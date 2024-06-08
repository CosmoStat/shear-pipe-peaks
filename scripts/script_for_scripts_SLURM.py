import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""                                                                
 
SCRIPT NAME:                                                                   
  script for generating slurm scripts for the data processing of cosmologies
  
outputs a bash script for each cosmology to be sent to the bash queue    
"""
)
parser.add_argument('job_name',type=str, help='job identifier')
parser.add_argument('work_dir',type=str, help='work directory for outputs')
parser.add_argument('run',type=str, help='run identifier')
parser.add_argument('n_cosmo',type=int, help='number of cosmologies')
parser.add_argument('n_tiles',type=int, help='number of tiles in footprint')
parser.add_argument('map_method',type=str, help='mass mapping method')
parser.add_argument('random_seed', type=int, nargs='?', default=None, help='random seed for noise generation. If not provided, random seed is None')
args = parser.parse_args()


job_name = args.job_name
work_dir = args.work_dir  
run = args.run                                                                 
n_cosmo = args.n_cosmo
map_method = args.map_method
n_tiles = args.n_tiles
random_seed = args.random_seed

# define directories
script_output_directory=work_dir+'/scripts_run_'+run+'_'+map_method+'/seed_'+str(random_seed)
# if the output directory is not present create it.
if not os.path.exists(script_output_directory): 
    os.makedirs(script_output_directory)

output_directory=work_dir+'/output_run_'+run+'_'+map_method+'/seed_'+str(random_seed)
# if the output directory is not present create it.
if not os.path.exists(output_directory): 
    os.makedirs(output_directory)


# make scripts
for cosmo in np.arange(n_cosmo):
    cosmo = str(cosmo).zfill(2)  # Converts cosmo to string and pads with leading zeros
    fileroot =  '%(script_output_directory)s/%(run)s%(cosmo)s'% locals()
    filename = fileroot+'.sh'
    print(filename)
    with open(filename, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --partition=comp\n")
        f.write("#SBATCH --output=%s.o\n" % fileroot)
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks-per-node=%d\n" % n_tiles)
        f.write("#SBATCH --time=5:00:00\n")
        f.write("#SBATCH --mem=32GB\n")
        f.write("#SBATCH --job-name=%s\n" % job_name)
        f.write("#SBATCH --mail-type=END\n")
        f.write("#SBATCH --mail-user=atersenov@physics.uoc.gr\n\n")
        
        f.write("echo 'This job runs on the following processors:'\n")
        f.write("NODES=$(scontrol show hostnames $SLURM_JOB_NODELIST)\n")
        f.write("echo $NODES\n")
        f.write("NPROCS=$(scontrol show nodes $SLURM_JOB_NODELIST | wc -l)\n")
        f.write("echo 'This job has allocated $NPROCS nodes'\n\n")

        f.write("module load healpix/3.82-ifx-2024.0\n\n")
        
        if random_seed is not None:
            f.write("~/miniconda3/envs/pysap/bin/python /home/tersenov/shear-pipe-peaks/scripts/data_processing_test.py %s /home/tersenov/shear-pipe-peaks/input/master_file.txt %s %d %s %s %s\n" % (cosmo, output_directory, n_tiles, map_method, run, random_seed))
        else:
            f.write("~/miniconda3/envs/pysap/bin/python /home/tersenov/shear-pipe-peaks/scripts/data_processing_test.py %s /home/tersenov/shear-pipe-peaks/input/master_file.txt %s %d %s %s\n" % (cosmo, output_directory, n_tiles, map_method, run))
        
        


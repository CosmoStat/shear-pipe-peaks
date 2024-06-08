import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""                                                                
 
SCRIPT NAME:                                                                   
  merged_cov_script_for_scripts.py

outputs a bash script for each noise realization to be sent to the bash queue    
"""
)
parser.add_argument('job_name',type=str, help='job identifier')
parser.add_argument('work_dir',type=str, help='work directory for outputs')
parser.add_argument('run',type=str, help='run identifier')
parser.add_argument('map_method',type=str, help='mass mapping method')
parser.add_argument('n_tiles',type=int, help='number of tiles in footprint')
parser.add_argument('n_realizations',type=int, help='number of noise realizations')
args = parser.parse_args()

# give meaningful variable names to the command line                           
# arguments:
job_name = args.job_name
work_dir = args.work_dir  
run = args.run                                                                 
map_method = args.map_method
n_tiles = args.n_tiles
n_realizations = args.n_realizations

# define directories
script_output_directory=work_dir+'/merged_scripts_run_'+run+'_'+map_method
output_directory=work_dir+'/merged_output_run_'+run+'_'+map_method
# if the output directory is not present create it.
if not os.path.exists(output_directory): 
    os.makedirs(script_output_directory)
    os.makedirs(output_directory)

# make scripts
for real in np.arange(n_realizations):
    real = str(real)
    fileroot =  '%(script_output_directory)s/%(run)s%(real)s'% locals()
    filename = fileroot+'.sh'
    print(filename)
    with open(filename, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --partition=pscomp\n")        
        f.write("#SBATCH --output=%s.o\n" % fileroot)
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks-per-node=%d\n" % n_tiles)
        f.write("#SBATCH --time=12:00:00\n")
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

        f.write("~/miniconda3/envs/pysap/bin/python /home/tersenov/shear-pipe-peaks/scripts/merged_bins_cov_cs.py %s /home/tersenov/shear-pipe-peaks/input/master_file_cov.txt %s %d %s %s\n" % (real, output_directory, n_tiles, map_method, run))
    
    


# ~/miniconda3/envs/pycs/bin/python /home/tersenov/shear-pipe-peaks/scripts/merged_bins_cov_cs.py %(real)s /home/tersenov/shear-pipe-peaks/input/master_file_cov.txt %(output_directory)s %(n_tiles)d %(map_method)s %(run)s  

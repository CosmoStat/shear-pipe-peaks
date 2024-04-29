import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""                                                                
 
SCRIPT NAME:                                                                   
  test.py

outputs a bash script for each cosmology to be sent to the bash queue    
"""
)
parser.add_argument('job_name',type=str, help='job identifier')
parser.add_argument('work_dir',type=str, help='work directory for outputs')
parser.add_argument('run',type=str, help='run identifier')
parser.add_argument('n_cosmo',type=int, help='number of cosmologies')
parser.add_argument('n_tiles',type=int, help='number of tiles in footprint')
parser.add_argument('map_method',type=str, help='mass mapping method')
args = parser.parse_args()

# give meaningful variable names to the command line                           
# arguments:
job_name = args.job_name
work_dir = args.work_dir  
run = args.run                                                                 
n_cosmo = args.n_cosmo
map_method = args.map_method
n_tiles = args.n_tiles

# define directories
script_output_directory=work_dir+'/scripts_run_'+run+'_'+map_method
output_directory=work_dir+'/output_run_'+run+'_'+map_method
# if the output directory is not present create it.
if not os.path.exists(output_directory): 
    os.makedirs(script_output_directory)
    os.makedirs(output_directory)

# make scripts
for cosmo in np.arange(n_cosmo):
    cosmo = str(cosmo).zfill(2)  # Converts cosmo to string and pads with leading zeros
    fileroot =  '%(script_output_directory)s/%(run)s%(cosmo)s'% locals()
    filename = fileroot+'.sh'
    print(filename)
    f = open(filename, 'w') 
    text = """#!/bin/bash                                                                    
#PBS -o %(fileroot)s.o 
#PBS -l nodes=1:ppn=%(n_tiles)d,walltime=5:00:00,mem=32GB
#PBS -N %(job_name)s                                                           
#PBS -j eo  
echo This jobs runs on the following processors:
NODES=`cat $PBS_NODEFILE`
echo $NODES
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes  

module load healpix/3.82-ifort-2023.0
module load gsl/2.7.1
module load fftw/3.3.9

~/miniconda3/envs/pycs/bin/python /home/tersenov/shear-pipe-peaks/scripts/merged_bins_data_processing.py %(cosmo)s /home/tersenov/shear-pipe-peaks/input/master_file.txt %(output_directory)s %(n_tiles)d %(map_method)s %(run)s  
    """ % locals()
    # write to file
    f.write(text)
    f.close()


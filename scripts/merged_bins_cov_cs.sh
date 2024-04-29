#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=19,walltime=13:00:00,mem=32GB
### job name
#PBS -N cov_paral_cs
### Redirect stdout and stderr to same file
#PBS -j oe
#PBS -m abe  
#PBS -M atersenov@physics.uoc.gr  

module load healpix/3.82-ifort-2023.0
module load gsl/2.7.1
module load fftw/3.3.9

## your bash script here:
~/miniconda3/envs/pycs/bin/python /home/tersenov/shear-pipe-peaks/scripts/merged_bins_cov_cs.py '/home/tersenov/shear-pipe-peaks' --nproc 19
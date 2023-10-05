#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=19,walltime=10:00:00,mem=32GB
### job name
#PBS -N l1-cov-par
### Redirect stdout and stderr to same file
#PBS -j oe
#PBS -m abe  
#PBS -M atersenov@physics.uoc.gr  

## your bash script here:
/home/tersenov/miniconda/bin/python /home/tersenov/shear-pipe-peaks/scripts/cov_l1_parallelized '/home/tersenov/shear-pipe-peaks' --nproc 19
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

## your bash script here:
/home/tersenov/.conda/envs/Sparse2D/bin/python /home/tersenov/shear-pipe-peaks/scripts/cov_paral_cs.py '/home/tersenov/shear-pipe-peaks' --nproc 19
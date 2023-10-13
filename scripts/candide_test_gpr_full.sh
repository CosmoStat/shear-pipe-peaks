#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=5,walltime=10:00:00,mem=32GB
### job name
#PBS -N test_gpr
### Redirect stdout and stderr to same file
#PBS -j oe
#PBS -m abe  
#PBS -M atersenov@physics.uoc.gr  

## your bash script here:
/home/tersenov/miniconda/envs/sppenv/bin/python /home/tersenov/shear-pipe-peaks/example/test_gpr_full.py '/home/tersenov/shear-pipe-peaks' --nproc 5

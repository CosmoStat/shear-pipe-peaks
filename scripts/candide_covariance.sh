#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=19,walltime=10:00:00,mem=32GB
### job name
#PBS -N peaks-covariance
### Redirect stdout and stderr to same file
#PBS -j oe

## your bash script here:
/home/tersenov/miniconda/bin/python /home/tersenov/shear-pipe-peaks/scripts/covariance.py '/home/tersenov/shear-pipe-peaks' --nproc 19
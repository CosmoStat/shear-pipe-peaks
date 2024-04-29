#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=48,walltime=10:00:00,mem=64GB
### job name
#PBS -N Inference
### Redirect stdout and stderr to same file
#PBS -j oe

## your bash script here:
~/miniconda3/envs/py36/bin/python /home/tersenov/shear-pipe-peaks/scripts/Inference.py '/home/tersenov/shear-pipe-peaks' --nproc 48
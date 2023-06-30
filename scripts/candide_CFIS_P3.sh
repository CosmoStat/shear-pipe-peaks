#!/bin/bash                                                                    
#PBS -k o
### resource allocation
#PBS -l nodes=1:ppn=8,walltime=10:00:00
### job name
#PBS -N peaks-mcmc
### Redirect stdout and stderr to same file
#PBS -j oe

## your bash script here:
~/.conda/envs/sp-peaks/bin/python /home/baumont/software/shear-pipe-peaks/example/constraints_CFIS-P3.py

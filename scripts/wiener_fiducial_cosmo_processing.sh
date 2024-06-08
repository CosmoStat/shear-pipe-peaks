#!/bin/bash
#SBATCH --partition=comp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=19
#SBATCH --time=8:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=wiener_fid
#SBATCH --mail-type=END
#SBATCH --mail-user=atersenov@physics.uoc.gr

echo 'This job runs on the following processors:'
NODES=$(scontrol show hostnames $SLURM_JOB_NODELIST)
echo $NODES
NPROCS=$(scontrol show nodes $SLURM_JOB_NODELIST | wc -l)
echo 'This job has allocated $NPROCS nodes'

module load healpix/3.82-ifx-2024.0

~/miniconda3/envs/pysap/bin/python /home/tersenov/shear-pipe-peaks/scripts/wiener_fiducial_cosmo_processing.py

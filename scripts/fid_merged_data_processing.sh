#!/bin/bash
#SBATCH --partition=pscomp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=19
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=mer_fid_cosmo_pr
#SBATCH --mail-type=END
#SBATCH --mail-user=atersenov@physics.uoc.gr

echo 'This job runs on the following processors:'
NODES=$(scontrol show hostnames $SLURM_JOB_NODELIST)
echo $NODES
NPROCS=$(scontrol show nodes $SLURM_JOB_NODELIST | wc -l)
echo 'This job has allocated $NPROCS nodes'

module load healpix/3.82-ifx-2024.0

~/miniconda3/envs/pysap/bin/python /home/tersenov/shear-pipe-peaks/scripts/fid_merged_data_processing.py

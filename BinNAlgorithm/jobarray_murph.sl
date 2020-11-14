#!/bin/bash
##SBATCH -N 1
#SBATCH --job-name=BinSims
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --array=0-1000
#SBATCH --output=./SLURMOUT/SLURM_log_%A-%a.out

module add r
R CMD BATCH "--no-save" mean_difference.R ./Rlogfiles/testingArray_$SLURM_ARRAY_TASK_ID.Rout

#!/bin/bash
##SBATCH -N 1
#SBATCH --job-name=MVNSTAN
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --time=01-00:00:00
#SBATCH --array=1-1000
#SBATCH --output=./SLURMOUT/slurm_log_%A-%a.out

module add r/3.6.0
module add gcc/6.3.0
R CMD BATCH "--no-save" CovarianceSimulations.R ./Rlogfiles/testingArray_$SLURM_ARRAY_TASK_ID.Rout

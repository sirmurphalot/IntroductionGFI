#!/bin/bash
##SBATCH -N 1
#SBATCH --job-name=RESim
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --time=05-00:00:00
#SBATCH --array=0-1000
#SBATCH --output=./SLURMOUT/murphmurphmurph_%A-%a.out

module add r/3.6.0
module add gcc/6.3.0
R CMD BATCH "--no-save" RESimulation.R ./Rlogfiles/testingArray_$SLURM_ARRAY_TASK_ID.Rout

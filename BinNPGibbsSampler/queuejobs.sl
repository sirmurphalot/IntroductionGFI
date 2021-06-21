#!/bin/bash
##SBATCH -N 1
#SBATCH --job-name=AllBinSim
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --time=06-00:00:00
#SBATCH --array=0-1199
#SBATCH --output=./SLURMOUT/YOURNAME_%A-%a.out

module add r/3.6.0
module add gcc/6.3.0
R CMD BATCH "--no-save" all_unknown_simulation.R ./Rlogfiles/testingArray_$SLURM_ARRAY_TASK_ID.Rout

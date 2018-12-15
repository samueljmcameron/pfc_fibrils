#!/bin/bash
#SBATCH --time=00:06:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out
#SBATCH --array=0-10

module restore standard_modules

python scan.py $SLURM_ARRAY_TASK_ID $1 $2

#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out
#SBATCH --array=0-7

module restore standard_modules


python scan_gamma.py $1 $2 $SLURM_ARRAY_TASK_ID

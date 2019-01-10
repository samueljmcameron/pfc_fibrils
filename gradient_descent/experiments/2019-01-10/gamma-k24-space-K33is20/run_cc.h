#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out
#SBATCH --array=0-100

module restore standard_modules


python scan2d.py $1 $2 $SLURM_ARRAY_TASK_ID

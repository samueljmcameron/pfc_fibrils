#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-arutenbe
#SBATCH --array=1,10-1000:10

module restore standard_modules


python pulling.py 0.12 1.0 $1 $SLURM_ARRAY_TASK_ID

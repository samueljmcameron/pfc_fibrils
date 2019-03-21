#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-arutenbe
#SBATCH --array=0-9

module restore standard_modules


python scanconstantgamma.py $1 $SLURM_ARRAY_TASK_ID 10.0

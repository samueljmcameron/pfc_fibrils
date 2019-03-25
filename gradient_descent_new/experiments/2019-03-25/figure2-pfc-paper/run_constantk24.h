#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-arutenbe
#SBATCH --array=0-4

module restore standard_modules


python scanconstantk24.py $SLURM_ARRAY_TASK_ID $1 10.0

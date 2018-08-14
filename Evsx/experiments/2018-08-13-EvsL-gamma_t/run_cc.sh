#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-arutenbe
#SBATCH --array 1-8

module restore standard_modules
module load python/2.7
module load scipy-stack

python generatedata_cc.py 1
#$SLURM_ARRAY_TASK_ID

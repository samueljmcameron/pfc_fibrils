#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%j.out

module restore standard_modules

python deltazero_energy.py $1 $2

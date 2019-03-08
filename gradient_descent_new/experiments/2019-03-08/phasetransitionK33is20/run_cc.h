#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=def-arutenbe

module restore standard_modules


python scan.py $1 $2 $3

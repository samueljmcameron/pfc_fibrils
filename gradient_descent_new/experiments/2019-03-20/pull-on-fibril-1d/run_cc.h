#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-arutenbe

module restore standard_modules


python pulling.py $1 $2 $3 $4

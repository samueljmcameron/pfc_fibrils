#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run-%j.out

module restore standard_modules

python raster_scan.py $1 $2

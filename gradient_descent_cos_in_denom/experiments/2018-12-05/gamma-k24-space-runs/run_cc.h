#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-arutenbe

module restore standard_modules

python raster_scan.py $1 $2

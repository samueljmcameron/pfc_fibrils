import numpy as np
import subprocess
import sys
from local_packages.singlerun import SingleRun

if __name__=="__main__":

    run = SingleRun("data/input.dat")

    run.run_exe()

    run.mv_file('psivsr')

    run.mv_file('xvals')

    run.concatenate_xvals()

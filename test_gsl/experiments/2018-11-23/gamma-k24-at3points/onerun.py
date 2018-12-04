import numpy as np
import subprocess
import sys
sys.path.append('../../gammak24_singlepoint/')
from singlerun import SingleRun

if __name__=="__main__":

    run = SingleRun("data/input.dat")

    run.run_exe()

    run.mv_file('psivsr')

    run.mv_file('observables')

    run.concatenate_observables()

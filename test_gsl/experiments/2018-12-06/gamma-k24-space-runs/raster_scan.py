import numpy as np
import subprocess
import sys
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun

if __name__=="__main__":


    datfile = "data/input.dat"
    executable = "../../../bin/gamma_k24_space"

    scan = {'\Lambda':sys.argv[1],'\omega':sys.argv[2]}
    run = SingleRun(datfile,executable=executable,scan=scan,
                    suffixlist=["K_{33}","\\Lambda","d_0","\\omega"])

    run.run_exe()

    run.mv_file('energy')
    run.mv_file('surfacetwist')
    run.mv_file('radius')
    run.mv_file('eta')
    run.mv_file('delta')


#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../bin')
from var_scan import one_run

if __name__=="__main__":

    # generates and then plots the data files for 
    # the specified parameter values

    eta_index = int(sys.argv[1])

    etas = np.logspace(-3,1,num=5)
    eta = etas[eta_index-1]

    one_run(eta)

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

    gamma_s_index = int(sys.argv[1])

    gamma_ss = np.logspace(-3,3,num=7)
    gamma_s = gamma_ss[gamma_s_index-1]

    one_run(gamma_s)

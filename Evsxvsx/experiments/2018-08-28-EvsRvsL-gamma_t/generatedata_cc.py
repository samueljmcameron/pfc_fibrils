#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../bin')
from var_scan2d import one_run


if __name__=="__main__":

    # generates and then plots the data files for 
    # the specified parameter values

    gamma_t_index = int(sys.argv[1])

    gamma_ts = np.logspace(-2,5,num=8)
    gamma_t = gamma_ts[gamma_t_index-1]

    one_run(gamma_t)

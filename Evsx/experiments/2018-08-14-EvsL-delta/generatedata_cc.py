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

    delta_index = int(sys.argv[1])

    deltas = np.logspace(-6,0,num=7)
    delta = deltas[delta_index-1]

    one_run(delta)

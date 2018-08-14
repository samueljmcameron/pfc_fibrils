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

    K33_index = int(sys.argv[1])

    K33s = np.logspace(-3,3,num=7)
    K33 = K33s[K33_index-1]

    one_run(K33)

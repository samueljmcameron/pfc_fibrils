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

    k24_index = int(sys.argv[1])

    k24s = np.linspace(-1,2,num=6,endpoint=True)
    k24 = k24s[k24_index-1]

    one_run(k24)

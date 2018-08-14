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

    Lambda_index = int(sys.argv[1])

    Lambdas = np.logspace(-1,6,num=8)
    Lambda = Lambdas[Lambda_index-1]

    one_run(Lambda)

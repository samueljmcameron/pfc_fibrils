#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../scripts')
from var_scan import one_run

if __name__=="__main__":

    # generates and then plots the data files for 
    # the specified parameter values

    #    Lambda_index = int(sys.argv[1])

    #    Lambdas = [0.01,0.1,1.0,10.0,20.0,30.0,40.0,50.0,60.0,
    #           70.0,80.0,90.0,100.0]#np.linspace(0,10,num=11,endpoint=True)
    #Lambda = Lambdas[Lambda_index-1]

    Lambda = float(sys.argv[1])
    one_run(Lambda)

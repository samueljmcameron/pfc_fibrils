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

    one_run(float(sys.argv[1]))

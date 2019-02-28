# testing the derivatives with this file turned out to be not so useful. see later tests.


import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
sys.path.append('../../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams
from energy import energy

from scipy.integrate import simps,romb


if __name__ == "__main__":

    if len(sys.argv)<5:

        user_input = input("input string of a gamma,k24,Lambda,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,Lambda,omega = user_input.split(',')

    else:
        gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]




        
    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\Lambda']=Lambda
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


    # read in file name info
    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)
        
    # create a class to do calculations with current parameters in scan.
    run = SingleRun(rp,executable="../../../bin/derivatives")
    # run C executable.
    run.run_exe()


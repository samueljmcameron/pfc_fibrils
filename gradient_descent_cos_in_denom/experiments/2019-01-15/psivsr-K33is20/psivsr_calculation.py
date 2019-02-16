################################################################################
# This file can only be run if the corresponding deltazero_energy.py script    #
# has been run first! Otherwise, 'observables_Emin' type file below will not   #
# exist and the script will fail. #
################################################################################

import numpy as np
import subprocess
import sys
import time
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":

    start_time = time.time()

    if len(sys.argv)<5:

        user_input = input("input string of a gamma,k24,Lambda,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,Lambda,omega = user_input.split(',')

    else:
        gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]


    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\Lambda']=str(Lambda)
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"]


    # read in file name info
    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)
        
    # create a class to do calculations with current parameters in scan.
    run = SingleRun(rp,executable="../../../bin/psivsr_calculation")

    # run C executable.
    run.run_exe()

    # move file written by C executable from temporary data path to true data path
    run.mv_file('observables')
    run.mv_file('psivsr')




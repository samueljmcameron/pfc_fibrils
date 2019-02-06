import numpy as np
import subprocess
import sys
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":


    gamma,k24 = sys.argv[1],sys.argv[2]

    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","\\omega","\\gamma_s"]
    
    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\omega']= '0'
    scan['\\Lambda']='0'



    # find minimum for delta = 0 case, so you know the upper bound for
    # the energy minimum.

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    run = SingleRun(rp)


    run.run_exe()

    run.mv_file('observables',newname='observables_Emin')
    

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
sys.path.append('../../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams
from energy import energy

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
    run = SingleRun(rp,executable="../../../bin/samplecalc_testEvsR_c")
    # run C executable.
    run.run_exe()
    dataname = run.mv_file('EvsR_c')

    data = np.loadtxt(dataname)


    R_cs = data[:,0]
    Es = data[:,1]
    dEdR_cs = data[:,2]

    i_max = np.argmax(Es)
    i_min = np.argmin(Es)

    fig,axarr = plt.subplots(2,sharex=True)

    axarr[0].plot(R_cs,Es,'--')
    axarr[0].scatter(R_cs[i_min],Es[i_min])
    axarr[0].scatter(R_cs[i_max],Es[i_max])
    axarr[1].plot(R_cs,dEdR_cs,'-.')
    axarr[1].set_ylim(-4,4)

    plt.show()



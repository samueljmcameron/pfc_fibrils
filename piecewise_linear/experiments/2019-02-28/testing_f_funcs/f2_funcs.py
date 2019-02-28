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
    run = SingleRun(rp,executable="../../../bin/f2_func")
    # run C executable.
    run.run_exe()
    dataname = run.mv_file('f2_func')

    data = np.loadtxt(dataname)


    zetas = data[:,0]
    f2s = data[:,1]
    df2dxis = data[:,2]
    df2dzetas = data[:,3]
        

    # plot Es vs psip_Rs and dEdpsip_Rs vs psip_Rs

    fig,axarr = plt.subplots(3,sharex=True)

    fig.set_size_inches(6,14)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    axarr[0].plot(zetas,f2s,'.',color=colors[0])
    axarr[0].set_ylabel(r"$f_1$")
    axarr[1].plot(zetas,df2dxis,'.',color=colors[1])
    axarr[1].set_ylabel(r"$\partial f_1/\partial \xi$")
    axarr[2].plot(zetas,df2dzetas,'.',color=colors[2])
    axarr[2].set_ylabel(r"$\partial f_1/\partial\zeta$")



    sname = dataname.replace("data/","results/").replace(".txt",".pdf")
    fig.subplots_adjust(left=0.2)

    plt.show()
    fig.savefig(sname)



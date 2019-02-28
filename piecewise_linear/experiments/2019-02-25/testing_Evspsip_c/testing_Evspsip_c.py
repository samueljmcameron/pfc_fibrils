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
    run = SingleRun(rp,executable="../../../bin/samplecalc_testEvspsip_c")
    # run C executable.
    run.run_exe()
    dataname = run.mv_file('Evspsip_c')

    data = np.loadtxt(dataname)


    psip_cs = data[:,0]
    Es = data[:,1]
    dEdpsip_cs = data[:,2]

    # determine where the min of Es occurs
    i_min = np.argmin(Es)

    # determine where the zero of dEdpsip_cs occur (should be at the same spot
    # as where the min above occurs according to calculus).
    j_zero1 = np.argmin(np.abs(dEdpsip_cs))
        
    dEs = np.gradient(Es,psip_cs)
    # plot Es vs psip_cs and dEdpsip_cs vs psip_cs

    fig,axarr = plt.subplots(2,sharex=True)

    axarr[0].plot(psip_cs,Es,'--')
    axarr[0].scatter(psip_cs[i_min],Es[i_min])
    axarr[0].set_ylabel(r"$E(psip_c)$")
    axarr[1].plot(psip_cs,dEdpsip_cs,'.')
    axarr[1].plot(psip_cs,dEs,'-')
    axarr[1].plot(psip_cs,np.abs(dEdpsip_cs),'--')
    axarr[1].scatter(psip_cs[j_zero1],dEdpsip_cs[j_zero1])
    axarr[1].set_ylabel(r"$dE/dpsip_c$")
    axarr[1].set_xlabel(r"$\psi_c^{\prime}$")



    print(f"minimum value of E = {Es[i_min]} is at psip_c = {psip_cs[i_min]}")

    # there should be a one to one correspondence between the above print
    # output and the below print output.

    print(f"first zero of dEdpsip_c = {dEdpsip_cs[j_zero1]} is at psip_c"
          f" = {psip_cs[j_zero1]}")

    sname = dataname.replace("data/","results/").replace(".txt",".pdf")

    plt.show()
    fig.savefig(sname)



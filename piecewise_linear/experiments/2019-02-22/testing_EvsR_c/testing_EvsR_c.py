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
    run = SingleRun(rp,executable="../../../bin/EvsR_c")
    # run C executable.
    run.run_exe()
    dataname = run.mv_file('EvsR_c')

    data = np.loadtxt(dataname)


    R_cs = data[:,0]
    Es = data[:,1]
    dEdR_cs = data[:,2]

    # determine where the max and mins of Es occur
    i_max = np.argmax(Es)
    i_min = np.argmin(Es)

    # determine where the zeros of dEdR_cs occur (should be at the same spot
    # as where the max and mins above occur according to calculus).
    j_zero1 = np.argmin(np.abs(dEdR_cs))
    
    j_zero2 = (np.argmin(np.abs(dEdR_cs[j_zero1+10:]))
               +j_zero1+10) # + 10 to ensure I search far enough away from the first minimum
    

    # plot Es vs R_cs and dEdR_cs vs R_cs

    fig,axarr = plt.subplots(2,sharex=True)

    axarr[0].plot(R_cs,Es,'--')
    axarr[0].scatter(R_cs[i_min],Es[i_min])
    axarr[0].scatter(R_cs[i_max],Es[i_max])
    axarr[0].set_ylabel(r"$E(R_c)$")
    axarr[1].plot(R_cs,dEdR_cs,'.')
    axarr[1].plot(R_cs,np.abs(dEdR_cs),'--')
    axarr[1].scatter(R_cs[j_zero1],dEdR_cs[j_zero1])
    axarr[1].scatter(R_cs[j_zero2],dEdR_cs[j_zero2])
    axarr[1].set_ylabel(r"$dE/dR_c$")
    axarr[1].set_xlabel(r"$R_c$")



    print(f"minimum value of E = {Es[i_min]} is at R_c = {R_cs[i_min]}")
    print(f"maximum value of E = {Es[i_max]} is at R_c = {R_cs[i_max]}")

    # there should be a one to one correspondence between the above two print
    # outputs and the below two print outputs.

    print(f"first zero of dEdR_c = {dEdR_cs[j_zero1]} is at R_c = {R_cs[j_zero1]}")
    print(f"second zero of dEdR_c = {dEdR_cs[j_zero2]} is at R_c = {R_cs[j_zero2]}") 

    sname = dataname.replace("data/","results/").replace(".txt",".pdf")

    fig.savefig(sname)
    plt.show()


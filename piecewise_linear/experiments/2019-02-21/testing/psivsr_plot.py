import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

from scipy.integrate import simps,romb

def psi(rs,R_c,R_s,psip_c,psip_s,psip_R):

    psi1 = (psip_c-psip_s)*R_c
    psi2 = (psip_s-psip_R)*R_s+psi1

    return np.piecewise(rs,[rs<=R_c,(rs<=R_s)&(rs>R_c),rs>R_s],
                        [lambda x: psip_c*x,lambda x: psip_s*x+psi1,
                         lambda x: psip_R*x+psi2])

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
    run = SingleRun(rp,executable="../../../bin/samplecalc")
    # run C executable.
    run.run_exe()

    # move file written by C executable from temporary data path to true data path
    run.mv_file('psivsr_discrete')



    c_data = np.loadtxt("data/_psivsr_3.0000e+01_1.0000e-01_1.0000e+03_5.0000e+00_4.0000e-02.txt")
    d_data = np.loadtxt("data/_psivsr_discrete_3.0000e+01_1.0000e-01_1.0000e+03_5.0000e+00_4.0000e-02.txt")



    c_rs = c_data[:,0]
    c_psis = c_data[:,1]
    c_rfs = c_data[:,3]

    c_dx = c_rs[1]-c_rs[0]

    d_rs = d_data[:,0]
    d_psis = d_data[:,1]
    d_rfs = d_data[:,3]

    d_dx = d_rs[1]-d_rs[0]

    R_c = float(run.params['R_c'])
    R_s = float(run.params['R_s'])
    R = float(run.params['R'])
    delta = float(run.params['delta'])
    gamma = float(run.params['\gamma_s'])
    omega = float(run.params['\omega'])
    k24 = float(run.params['k_{24}'])

    
    E = 2/(R*R)*simps(d_rfs,d_rs)

    E += omega*delta*delta/2*(3/4*delta*delta-1)
    E += -(1+k24)/(R*R)*np.sin(c_psis[-1])+2*gamma/R

    print("E python discrete from continuum integral = ",E)

    """
    fig, axarr = plt.subplots(2,sharex=True)
    
    axarr[0].plot(c_rs,c_psis,'o',label='continuous')
    axarr[0].plot(d_rs,d_psis,'.',label='discrete')
    axarr[0].set_ylabel(r"$\psi(r)$")

    axarr[1].plot(c_rs,c_rfs,'o',label='continuous')
    axarr[1].plot(d_rs,d_rfs,'.',label='discrete')
    axarr[1].set_ylabel(r"$rf(r)$")
    axarr[1].set_xlabel(r"$r$")

    plt.show()
    """

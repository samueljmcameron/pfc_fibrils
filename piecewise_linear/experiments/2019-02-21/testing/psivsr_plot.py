import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

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



    R_c = 0.0355
    R_s = 0.865
    psip_c = 3
    psip_s = 0.0
    psip_R = 0.5

        
    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\Lambda']=Lambda
    scan['\\omega']= omega
    scan['R_c']=str(R_c)
    scan['R_s']=str(R_s)
    scan['psip_c']=str(psip_c)
    scan['psip_s']=str(psip_s)
    scan['psip_R']=str(psip_R)
    
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

    d_rs = d_data[:,0]
    d_psis = d_data[:,1]

    plt.plot(c_rs,c_psis,'o')

    plt.plot(d_rs,d_psis,'.')
    plt.plot(c_rs,psi(c_rs,R_c,R_s,psip_c,psip_s,psip_R),'k-')

    plt.show()

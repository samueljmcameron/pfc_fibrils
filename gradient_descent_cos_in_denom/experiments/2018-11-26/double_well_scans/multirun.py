################################################################################
# In this file, generate E vs t files for multiple values of omega, which are  #
# read in from .dat files in data/ folder. #
################################################################################




import numpy as np
import sys

sys.path.append('../../local_packages/')
from doublewell import DoubleWell

def read_in_xs(data):
    omegas = data[:,0]
    Lambdas = data[:,1]
    Es = data[:,2]
    Rs = data[:,3]
    etas = data[:,4]
    deltas = data[:,5]

    return omegas,Lambdas,Es,Rs,etas,deltas

def write_dict(s,omega,Lambda,R0,R1,eta0,eta1,delta0,delta1):

    s['\omega'] = str(omega)
    s['\Lambda'] = str(Lambda)
    s['R0'] = str(R0)
    s['R1'] = str(R1)
    s['eta0'] = str(eta0)
    s['eta1'] = str(eta1)
    s['delta0'] = str(delta0)
    s['delta1'] = str(delta1)

    return



if __name__=="__main__":

    front_scans = np.loadtxt("data/front_scans.dat")
    
    back_scans = np.loadtxt("data/back_scans.dat")

    omegas,Lambdas,Es,R0s,eta0s,delta0s=read_in_xs(front_scans)
        
    omegas,Lambdas,Es,R1s,eta1s,delta1s=read_in_xs(back_scans)

    s = {}

    for i,omega in enumerate(omegas):

        write_dict(s,omegas[i],Lambdas[i],R0s[i],R1s[i],eta0s[i],
                   eta1s[i],delta0s[i],delta1s[i])
    
        dw = DoubleWell(scan=s);

        dw.run_exe()

        dw.mv_file()

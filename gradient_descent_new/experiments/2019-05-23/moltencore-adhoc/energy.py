import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import sys
sys.path.append('../../scripts/')
from psidata import PsiData
from observabledata import ObservableData



def rf_frank(r,psi,psiprime,K33):

    a = np.where(r != 0.0,0.5*(psiprime+0.5*np.sin(2*psi)/r-1)**2
                 +0.5*K33*np.sin(psi)**4/r**2,0)

    return r*a


def rf_pf(r,psi,psiprime,delta,eta,Lambda):

    a = Lambda*delta**2/4.0*(4*np.pi*np.pi-eta*eta*np.cos(psi)**2)**2

    return r*a

def rf_total(r,psi,psiprime,delta,eta,K33,Lambda):

    return rf_frank(r,psi,psiprime,K33)+rf_pf(r,psi,psiprime,delta,eta,Lambda)


def E(R0,rs,psis,psiprimes,delta,eta,K33,k24,Lambda,omega,gamma):

    i_R0 = len(rs[rs<R0])

    R = rs[-1]

    E = 2.0/R**2*simps(rf_frank(rs,psis,psiprimes,K33),rs)

    E += 2.0/R**2*simps(rf_pf(rs[i_R0:],psis[i_R0:],psiprimes[i_R0:],
                              delta,eta,Lambda),rs[i_R0:])


    E += delta**2*omega*(0.75*delta*delta-1)*0.5

    E += 1/R*(-(1+k24)*np.sin(psis[-1])**2/R+2.0*gamma)

    return E

if __name__ == "__main__":

    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf=loadsuf


    gamma = 0.04
    k24 = 0.5
    Lambda = 600
    omega = 20
    K33 = 30.0


    scan = {}
    scan['\\gamma_s'] = str(gamma)
    scan['k_{24}'] = str(k24)
    scan['\\Lambda'] = str(Lambda)
    scan['\\omega'] = str(omega)
    scan['K_{33}'] = str(K33)

    observablestuff = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                                     name=f"observables_frustrated")

    delta = observablestuff.delta()
    eta = observablestuff.eta()
    E0 = observablestuff.E()

    print(delta,eta,E0)


    psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,name=f"psivsr_frustrated",
                       sfile_format="pdf")

    rs = psistuff.r()
    psis = psistuff.psi()
    psiprimes = psistuff.psiprime()

    R0s = np.copy(rs[:len(rs)-10:5])

    Es = np.copy(R0s)*0


    for i,R0 in enumerate(R0s):

        Es[i] = E(R0,rs,psis,psiprimes,delta,eta,K33,k24,Lambda,omega,gamma)


    plt.plot(R0s,Es,'.')

    plt.plot(R0s,E0+R0s*0,'--')

    plt.show()

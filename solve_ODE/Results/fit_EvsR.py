from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
from scipy.optimize import curve_fit

font = {'family':'sans-serif','weight':'normal','size':'22'}
rc('font', **font)
rc('text',usetex=True)
params = {'text.latex.preamble':[r'\usepackage{siunitx}',
                                 r'\usepackage{sfmath}',
#                                 r'\sisetup{detect-family}=true',
                                 r'\usepackage{amsmath}']}

##########################################################################################
# Make a plot of E(qR) and \psi(qR). Note: R != R_{\text{min}}.                          #
##########################################################################################


def p_R(R,k24,K33):
    return 2*(1-k24)/(K33*R*R)

def q_R(R,K33):
    return 2.0/(K33*R*R)

def psi_prime_lin(R,k24,K33):
    # assuming psi(r)~psi'_0*r, determine psi'_0
    # which minimizes the free energy in the small
    # angle approximation
    p = 2*(1-k24)/(K33*R*R)
    q = 2.0/(K33*R*R)
    w = np.cbrt(0.5*(q+np.sqrt(q*q+4.0/27.0*p*p*p)))
    return w-p/(3.0*w)

def E_scaled_lin(R,k24,K33,gamma):
    # assuming psi(r)~psi'_0*r, calculate the integrated free
    # energy per unit volume of the fibril phase
    pp = psi_prime_lin(R,k24,K33)
    return 0.5*(1-2*pp)**2-(k24+1)*pp**2+K33/4.0*pp**4*R**2+2*gamma/R

def func(qR,invariant,qReq,alpha):
    # fitted functional form for E
    sigma = 2*invariant/(alpha*(qReq)**(1-alpha))
    return 2*invariant/qR-sigma/qR**(alpha)+1.0/2.0

def find_qReq_alpha(qR,EK22q2,invariant):
    # determine qReq parameter from data by simply 
    # choosing qR[i] such that EK22q2[i] is the min
    # value of EK22q2
    qReq = qR[np.argmin(EK22q2)]
    # use curve_fit to determine alpha parameter
    popt, pcov = curve_fit(lambda qR,
                           alpha:func(qR,invariant,qReq,alpha),
                           qR,EK22q2)
    alpha = popt[0]
    return qReq,alpha


plt.rcParams.update(params)

K33_K22=30.0
K24_K22 = float(sys.argv[1])
invariant = float(sys.argv[2])
qoverPI = float(sys.argv[3])

k24 = 2*K24_K22-1


savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/EvsR.c

fbase="qoverPIis%1.4e_K24K22is_%1.4e_INVARIANTis%1.4eEvsR"%(qoverPI,K24_K22,invariant)

# first load and plot E vs R, and psi(R) vs R

EvsR = np.loadtxt(fbase + ".txt")
qR = EvsR[:,0]
EK22q2 = EvsR[:,1]

# find qReq, alpha

qReq,alpha = find_qReq_alpha(qR,EK22q2,invariant)

fig1,ax = plt.subplots()
fig1.set_size_inches(12,12)
ax.plot(qR,EK22q2,'.',label = r"$%1.3lf,%1.3lf$"\
                   %(K24_K22,invariant))
ax.plot(qR,E_scaled_lin(qR,k24,K33_K22,invariant),'-')
ax.plot(qR,func(qR,invariant,qReq,alpha),
        '--',label = r"$\alpha=%0.2f,\:qR_{\text{eq}}=%.3f$"%(alpha,qReq))
ax.set_ylabel(r"$E/(K_{22}q^2)$")
ax.set_xlabel(r"$qR$")
ax.set_xscale(u'log')
#ax.set_yscale(u'log')
legend = ax.legend(loc='best',title = "\t" r"$\frac{K_{24}}{K_{22}}$" ", " 
                          r"$\frac{\gamma}{K_{22}q}$")
legend.get_title().set_fontsize(20)


fig1.tight_layout()


fig1.savefig(savePath + "EvsRfits/" + fbase + ".pdf",dpi = 'figure')

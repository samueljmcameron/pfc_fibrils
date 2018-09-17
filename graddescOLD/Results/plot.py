from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc

font = {'family':'sans-serif','weight':'normal','size':'22'}
rc('font', **font)
rc('text',usetex=True)
params = {'text.latex.preamble':[r'\usepackage{siunitx}',
                                 r'\usepackage{sfmath}',
#                                 r'\sisetup{detect-family}=true',
                                 r'\usepackage{amsmath}']}

##########################################################################################
# Make a plot of E(R) and \psi(R). Note: R != R_{\text{min}}.                          #
##########################################################################################


plt.rcParams.update(params)


k24 = float(sys.argv[1])
gamma = float(sys.argv[2])
Lambda = float(sys.argv[3])
plotE = sys.argv[4]
plotPsi = sys.argv[5]

savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/EvsR.c

fbase="k24_%1.4e_gamma_s%1.4e_Lambda_%1.4e_EvsR"%(k24,gamma,Lambda)
Psibase="k24_%1.4e_gamma_s%1.4e_Lambda_%1.4e_PsiVSr"%(k24,gamma,Lambda)

# first load and plot E vs R, and psi(R) vs R

if plotE == "yes":
    
    EvsR = np.loadtxt(fbase + ".txt")
    R = EvsR[:,0]
    E = EvsR[:,1]
    dEdR = EvsR[:,2]
    psiR = EvsR[:,3]
    
    
    fig1,axarr1 = plt.subplots(3,sharex = True)
    fig1.set_size_inches(8,12)
    axarr1[0].plot(R,E,'-',label = r"$%1.3lf,%1.3lf$"\
                       %(k24,gamma))
    axarr1[0].set_ylabel(r"$E$")
    axarr1[0].set_xlabel(r"$R$")
    #    axarr1[0].set_xscale(u'log')
    legend = axarr1[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    legend.get_title().set_fontsize(20)
    axarr1[1].plot(R,dEdR,'-')
    axarr1[1].set_ylabel(r"$\frac{dE}{dR}$")

    axarr1[2].plot(R,psiR,'-')
    #axarr1[2].plot(R,psiRapprox(k24,gamma,R),'--')
    axarr1[2].set_ylabel(r"$\psi^{\prime}(R)$")
    axarr1[2].set_xlabel(r"$R$")
    #    axarr1[2].set_xscale(u'log')
    
    fig1.tight_layout()
    fig1.savefig(savePath + "EvsRgraphs/" + fbase + ".pdf")

##########################################################################################
# Make a plot of \psi(r) and \psi'(R).                                                 #
##########################################################################################

def psi_ps(r,dpsi0,K33):
    return dpsi0*r+((3*K33-4)*dpsi0**3-3*dpsi0**2)/12*r**3

def dpsi_ps(r,dpsi0,K33):
    return dpsi0+3*((3*K33-4)*dpsi0**3-3*dpsi0**2)/12*r**2

if plotPsi == "yes":
    
    PsiVSr = np.loadtxt(Psibase + ".txt")
    r = PsiVSr[:,0]
    psi = PsiVSr[:,1]
    dpsi_dr_q = PsiVSr[:,2]
    
    fig2,axarr2 = plt.subplots(2,sharex = True)
    fig2.set_size_inches(12,12)
    axarr2[0].plot(r,psi,'-',label = r"$%1.3lf,%1.3lf$"\
                       %(k24,gamma))
#    axarr2[0].plot(r,psi_ps(r,dpsi_dr_q[0],30),'--')
    axarr2[0].set_ylabel(r"$\psi(r)$" " " r"$(\si{\radian})$")
    axarr2[0].set_xlabel(r"$r$")
    legend = axarr2[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    legend.get_title().set_fontsize(20)
    axarr2[1].plot(r,dpsi_dr_q,'-')
#    axarr2[1].plot(r,dpsi_ps(r,dpsi_dr_q[0],30),'--')
    axarr2[1].set_ylabel(r"$\psi^{\prime}(r)$")
    axarr2[1].set_xlabel(r"$r$")
    fig2.savefig(savePath + "PsiGraphs/" + Psibase + ".pdf")

plt.show()

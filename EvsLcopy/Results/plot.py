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
# Make a plot of E(L) and \psi(R). Note: R != R_{\text{min}}.                          #
##########################################################################################


plt.rcParams.update(params)


d0 = float(sys.argv[1])
R = float(sys.argv[2])
eta = float(sys.argv[3])
plotE = sys.argv[4]
plotPsi = sys.argv[5]

savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/main.c


fbase="d0_%1.4e_R_%1.4e_eta_%1.4e_EvsL",%(d0,R,eta)
Psibase="d0_%1.4e_R_%1.4e_eta_%1.4e_PsiVSr",%(d0,R,eta)

# first load and plot E vs L, and psi(R) vs L

if plotE == "yes":
    
    EvsL = np.loadtxt(fbase + ".txt")
    L = EvsL[:,0]
    E = EvsL[:,1]
    dEdL = EvsL[:,2]
    psiL = EvsL[:,3]
    
    
    fig1,axarr1 = plt.subplots(3,sharex = True)
    fig1.set_size_inches(8,12)
    axarr1[0].plot(L,E,'-',label = r"$%1.3lf,%1.3lf$"\
                       %(k24,gamma))
    axarr1[0].set_ylabel(r"$E$")
    axarr1[0].set_xlabel(r"$L$")
    #    axarr1[0].set_xscale(u'log')
    legend = axarr1[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    legend.get_title().set_fontsize(20)
    axarr1[1].plot(L,dEdL,'-')
    axarr1[1].set_ylabel(r"$\frac{dE}{dL}$")

    axarr1[2].plot(L,psiL,'-')
    axarr1[2].set_ylabel(r"$\psi(R)$")
    axarr1[2].set_xlabel(r"$L$")
    #    axarr1[2].set_xscale(u'log')
    
    fig1.tight_layout()
    fig1.savefig(savePath + "EvsLgraphs/" + fbase + ".pdf")

##########################################################################################
# Make a plot of \psi(r) and \psi'(R).                                                 #
##########################################################################################

if plotPsi == "yes":
    
    PsiVSr = np.loadtxt(Psibase + ".txt")
    r = PsiVSr[:,0]
    psi = PsiVSr[:,1]
    dpsi_dr = PsiVSr[:,2]
    
    fig2,axarr2 = plt.subplots(2,sharex = True)
    fig2.set_size_inches(12,12)

    length = 2*2*2*2*2*2*2*2*2*2*2+1
    i = 1
    while (len(r)-i*length>0):        
        axarr2[0].plot(r[(i-1)*length:i*length],psi[(i-1)*length:i*length],
                       '-',label = r"$%d$"\
                       %(i))
        axarr2[1].plot(r[(i-1)*length:i*length],dpsi_dr[(i-1)*length:i*length],'-')
        i += 1
    axarr2[0].set_ylabel(r"$\psi(r)$" " " r"$(\si{\radian})$")
    axarr2[0].set_xlabel(r"$r$")
    legend = axarr2[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    legend.get_title().set_fontsize(20)
    axarr2[1].set_ylabel(r"$\psi^{\prime}(r)$")
    axarr2[1].set_xlabel(r"$r$")
    fig2.savefig(savePath + "PsiGraphs/" + Psibase + ".pdf")

plt.show()

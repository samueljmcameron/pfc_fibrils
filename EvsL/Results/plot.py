from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../bin')
from fig_settings import configure_fig_settings
import seaborn.apionly as sns


configure_fig_settings()
colors = sns.color_palette('muted')

##########################################################################################
# Make a plot of E(R) and \psi(R). Note: R != R_{\text{min}}.                          #
##########################################################################################

k24 = 0.8
gamma = 0.05


d0 = float(sys.argv[1])
L = float(sys.argv[2])
eta = float(sys.argv[3])
plotE = sys.argv[4]
plotpsi = sys.argv[5]

savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/main.c

datapath = "../data/"
fbase="d0_%1.4e_L_%1.4e_eta_%1.4e_EvsR"%(d0,L,eta)
psibase="d0_%1.4e_L_%1.4e_eta_%1.4e_psiVSr_R"%(d0,L,eta)

# first load and plot E vs R, and psi(R) vs R

if plotE == "yes":
    
    EvsR = np.loadtxt(datapath + fbase + ".txt")
    R = EvsR[:,0]
    E = EvsR[:,1]
    dEdR = EvsR[:,2]
    psiR = EvsR[:,3]
    
    
    fig1,axarr1 = plt.subplots(3,sharex = True)
    fig1.set_size_inches(8,12)
    axarr1[0].plot(R,E,'.',color = colors[0])
    axarr1[0].set_ylabel(r"$E$")
    axarr1[0].set_xlabel(r"$R$")
    #    axarr1[0].set_xscale(u'log')
    legend = axarr1[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    axarr1[1].plot(R,dEdR,'.',color = colors[1])
    axarr1[1].set_ylabel(r"$\frac{dE}{dR}$")

    axarr1[2].plot(R,psiR,'.',color = colors[2])
    axarr1[2].set_ylabel(r"$\psi(R)$")
    axarr1[2].set_xlabel(r"$R$")
    #    axarr1[2].set_xscale(u'log')
    
    fig1.tight_layout()
    fig1.savefig(savePath + fbase + ".pdf")

##########################################################################################
# Make a plot of \psi(r) and \psi'(R).                                                 #
##########################################################################################

if plotpsi == "yes":
    
    psiVSr = np.loadtxt(psibase + ".txt")
    r = psiVSr[:,0]
    psi = psiVSr[:,1]
    dpsi_dr = psiVSr[:,2]
    
    fig2,axarr2 = plt.subplots(2,sharex = True)
    fig2.set_size_inches(12,12)

    length = 2*2*2*2*2*2*2*2*2*2*2+1
    i = 1
    while (len(r)-i*length>0):        
        axarr2[0].plot(r[(i-1)*length:i*length],psi[(i-1)*length:i*length],
                       '.',color = colors[i%len(colors)],label = r"$%d$"\
                           %(i))
        axarr2[1].plot(r[(i-1)*leength:i*length],dpsi_dr[(i-1)*length:i*length],
                       '.',color = colors[i%len(colors)])
        i += 1
    axarr2[0].set_ylabel(r"$\psi(r)$" " " r"$(\si{\radian})$")
    axarr2[0].set_xlabel(r"$r$")
    legend = axarr2[0].legend(loc='best',title = "\t" r"$k_{24}$" ", " 
                              r"$\gamma$")
    legend.get_title().set_fontsize(20)
    axarr2[1].set_ylabel(r"$\psi^{\prime}(r)$")
    axarr2[1].set_xlabel(r"$r$")
    fig2.savefig(savePath + psibase + ".pdf")

plt.show()

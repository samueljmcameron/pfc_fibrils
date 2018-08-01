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
# Make a plot of E(Rq) and \psi(Rq). Note: R != R_{\text{min}}.                          #
##########################################################################################


plt.rcParams.update(params)


k24 = 0.8
Lambda = float(sys.argv[1])
R = float(sys.argv[2])
plotE = sys.argv[3]
plotpsi = sys.argv[4]
plotrf = sys.argv[5]

savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/EvsR.c

fbase="E_vs_it_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R)
psibase="psi_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R)
rfbase="rf_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R)

# first load and plot E vs R, and psi(R) vs R

if plotE == "yes":
    EVSit = np.loadtxt(fbase + ".txt")
    it = EVSit[:,0]
    E = EVSit[:,1]
    

    fig1,ax1 = plt.subplots()
    fig1.set_size_inches(12,12)
    ax1.plot(it,E,'.',label = r"$%1.3lf,%1.3lf$"\
                 %(k24,Lambda))
    ax1.plot(it,E*0 + E[-1],'--')
    ax1.set_ylabel(r"$E$")
    ax1.set_xlabel(r"$t$")

    fig1.tight_layout()
    fig1.savefig(savePath + "EVSit/" + fbase + ".pdf")

##########################################################################################
# Make a plot of \psi(rq) and \psi'(Rq).                                                 #
##########################################################################################

if plotpsi == "yes":
    
    PsiVSr = np.loadtxt(psibase + ".txt")
    r = PsiVSr[:,0]
    psi = PsiVSr[:,1]
    dpsi_dr = PsiVSr[:,2]
    
    fig2,axarr2 = plt.subplots(2,sharex=True)
    fig2.set_size_inches(12,12)
    axarr2[0].plot(r,psi,'.',label = r"$%1.3lf,%1.3lf$"\
                       %(k24,Lambda))
    axarr2[0].set_ylabel(r"$\psi(r)$")
    axarr2[0].set_yscale('log')
    axarr2[1].plot(r,dpsi_dr,'.',label = r"$%1.3lf,%1.3lf$"\
                       %(k24,Lambda))
    axarr2[1].set_ylabel(r"$\psi'(r)$")
    axarr2[1].set_xlabel(r"$r$")
    axarr2[1].set_yscale('log')
    fig2.savefig(savePath + "psiVSr/" + psibase + ".pdf")

if plotrf == "yes":
    rfVSr = np.loadtxt(rfbase + ".txt")
    r = rfVSr[:,0]
    rf = rfVSr[:,1]
    

    fig1,ax1 = plt.subplots()
    fig1.set_size_inches(12,12)
    ax1.plot(r,rf,'.',label = r"$%1.3lf,%1.3lf$"\
                 %(k24,Lambda))
    ax1.set_ylabel(r"$rf$")
    ax1.set_xlabel(r"$r$")
    ax1.set_yscale('log')
    fig1.tight_layout()
    fig1.savefig(savePath + "rfVSr/" + rfbase + ".pdf")



plt.show()

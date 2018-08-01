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
R1 = float(sys.argv[2])
R2 = float(sys.argv[3])
plotE = sys.argv[4]
plotpsi = sys.argv[5]
plotrf = sys.argv[6]

savePath="Figures/" # figures will be saved in subfolder

# define all the saving/loading text stuff from ../code/EvsR.c

fbase1="E_vs_it_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R1)
psibase1="psi_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R1)
rfbase1="rf_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R1)

fbase2="E_vs_it_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R2)
psibase2="psi_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R2)
rfbase2="rf_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e"%(k24,Lambda,R2)

# first load and plot E vs R, and psi(R) vs R

if plotE == "yes":
    EVSit1 = np.loadtxt(fbase1 + ".txt")
    it1 = EVSit1[:,0]
    E1 = EVSit1[:,1]
    EVSit2 = np.loadtxt(fbase2 + ".txt")
    it2 = EVSit2[:,0]
    E2 = EVSit2[:,1]
    

    fig1,ax1 = plt.subplots()
    fig1.set_size_inches(12,12)
    ax1.plot(it1,E1,'b.',label = r"$%1.3lf,%1.3lf$"
             %(k24,Lambda))
    ax1.plot(it2,E2,'r.',label = r"$%1.3lf,%1.3lf$"
             %(k24,Lambda))
    ax1.set_ylabel(r"$E$")
    ax1.set_xlabel(r"$t$")
#    ax1.set_yscale('log')

    fig1.tight_layout()
    fig1.savefig(savePath + "compare/" + fbase1 + ".pdf")

##########################################################################################
# Make a plot of \psi(rq) and \psi'(Rq).                                                 #
##########################################################################################

if plotpsi == "yes":
    
    PsiVSr1 = np.loadtxt(psibase1 + ".txt")
    r1 = PsiVSr1[:,0]
    psi1 = PsiVSr1[:,1]
    dpsi_dr1 = PsiVSr1[:,2]

    PsiVSr2 = np.loadtxt(psibase2 + ".txt")
    r2 = PsiVSr2[:,0]
    psi2 = PsiVSr2[:,1]
    dpsi_dr2 = PsiVSr2[:,2]
    
    fig2,axarr2 = plt.subplots(2,sharex=True)
    fig2.set_size_inches(12,12)
    axarr2[0].plot(r1,psi1,'b.',label = r"$%1.3lf,%1.3lf$"
                   %(k24,Lambda))
    axarr2[0].plot(r2,psi2,'r.',label = r"$%1.3lf,%1.3lf$"
                   %(k24,Lambda))
    axarr2[0].set_ylabel(r"$\psi(r)$")
    #    axarr2[0].set_yscale('log')
    axarr2[1].plot(r1,dpsi_dr1,'b.',label = r"$%1.3lf,%1.3lf$"
                   %(k24,Lambda))
    axarr2[1].plot(r2,dpsi_dr2,'r.',label = r"$%1.3lf,%1.3lf$"
                   %(k24,Lambda))
    axarr2[1].set_ylabel(r"$\psi'(r)$")
    axarr2[1].set_xlabel(r"$r$")
    #    axarr2[1].set_yscale('log')
    fig2.savefig(savePath + "compare/" + psibase1 + ".pdf")

if plotrf == "yes":
    rfVSr1 = np.loadtxt(rfbase1 + ".txt")
    r1 = rfVSr1[:,0]
    rf1 = rfVSr1[:,1]
    rfVSr2 = np.loadtxt(rfbase2 + ".txt")
    r2 = rfVSr2[:,0]
    rf2 = rfVSr2[:,1]
    

    fig1,ax1 = plt.subplots()
    fig1.set_size_inches(12,12)
    ax1.plot(r1,rf1,'b.',label = r"$%1.3lf,%1.3lf$"
             %(k24,Lambda))
    ax1.plot(r2,rf2,'r.',label = r"$%1.3lf,%1.3lf$"
             %(k24,Lambda))
    ax1.set_ylabel(r"$rf_1$")
    ax1.set_xlabel(r"$r$")
    #    ax1.set_yscale('log')
    fig1.tight_layout()
    fig1.savefig(savePath + "compare/" + rfbase1 + ".pdf")



plt.show()

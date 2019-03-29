import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


data1 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+00_2.5000e+01_4.0000e-02.txt")
data2 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+01_2.5000e+01_4.0000e-02.txt")
data3 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+02_2.5000e+01_4.0000e-02.txt")
data4 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+03_2.5000e+01_4.0000e-02.txt")


strains1 = data1[:,0]*100
energys1 = data1[:,1]
stresses1 = np.gradient(energys1,strains1/100.0)
stresses1[0] = 0.0
Rs1 = data1[:,2]
deltas1 = data1[:,4]
twists1 = data1[:,5]


strains2 = data2[:,0]*100
energys2 = data2[:,1]
stresses2 = np.gradient(energys2,strains2/100.0)
stresses2[0] = 0.0
Rs2 = data2[:,2]
deltas2 = data2[:,4]
twists2 = data2[:,5]

strains3 = data3[:,0]*100
energys3 = data3[:,1]
stresses3 = np.gradient(energys3,strains3/100.0)
stresses3[0] = 0.0
Rs3 = data3[:,2]
deltas3 = data3[:,4]
twists3 = data3[:,5]

strains4 = data4[:,0]*100
energys4 = data4[:,1]
stresses4 = np.gradient(energys4,strains4/100.0)
stresses4[0] = 0.0
Rs4 = data4[:,2]
deltas4 = data4[:,4]
twists4 = data4[:,5]



fig0,axarr0 = plt.subplots(4,sharex=True)
fig0.set_size_inches(5,4*4)

axarr0[0].plot(strains1,energys1,'-',label=r"$\Lambda=1$")
axarr0[0].plot(strains2,energys2,'-.',label=r"$\Lambda=10$")
axarr0[0].plot(strains3,energys3,'.-',label=r"$\Lambda=100$")
axarr0[0].plot(strains4,energys4,'.-',label=r"$\Lambda=1000$")

axarr0[0].set_ylabel(r"$E$")

axarr0[1].plot(strains1,deltas1/deltas1[0],'-',label=r"$\Lambda=1$")
axarr0[1].plot(strains2,deltas2/deltas2[0],'-.',label=r"$\Lambda=10$")
axarr0[1].plot(strains3,deltas3/deltas3[0],'.-',label=r"$\Lambda=100$")
axarr0[1].plot(strains4,deltas4/deltas4[0],'.-',label=r"$\Lambda=1000$")

axarr0[1].set_ylabel(r"$\delta/\delta_0$")

axarr0[2].plot(strains1,twists1,'-',label=r"$\Lambda=1$")
axarr0[2].plot(strains2,twists2,'-.',label=r"$\Lambda=10$")
axarr0[2].plot(strains3,twists3,'.-',label=r"$\Lambda=100$")
axarr0[2].plot(strains4,twists4,'.-',label=r"$\Lambda=1000$")

axarr0[2].set_ylabel(r"$\psi(R)$")

axarr0[3].plot(strains1,Rs1,'-',label=r"$\Lambda=1$")
axarr0[3].plot(strains2,Rs2,'-.',label=r"$\Lambda=10$")
axarr0[3].plot(strains3,Rs3,'.-',label=r"$\Lambda=100$")
axarr0[3].plot(strains4,Rs4,'.-',label=r"$\Lambda=1000$")

axarr0[3].set_ylabel(r"$R$" + " (poisson ratio of 0.5)")


axarr0[3].set_xlabel("strain (%)")

axarr0[1].legend(frameon=False,loc="lower left")


fig0.subplots_adjust(left=0.2,right=0.95)
fig0.savefig("results/roughplot.pdf")

fig1,ax1 = plt.subplots()
fig1.set_size_inches(5,4)


ax1.plot(strains1,stresses1,'-',label=r"$\Lambda=1$")
ax1.plot(strains2,stresses2,'-.',label=r"$\Lambda=10$")
ax1.plot(strains3,stresses3,'.-',label=r"$\Lambda=100$")
ax1.plot(strains4,stresses4,'.-',label=r"$\Lambda=1000$")
ax1.set_ylabel(r"$\sigma$")
ax1.set_xlabel("strain (%)")
ax1.set_ylim(bottom=0)
ax1.legend(frameon=False,loc="upper left")
axins = inset_axes(ax1,width=1.3,height=0.9)

axins.plot(strains1[:10],stresses1[:10],'-',label=r"$\Lambda=1$")
axins.plot(strains2[:10],stresses2[:10],'-.',label=r"$\Lambda=10$")
axins.plot(strains3[:10],stresses3[:10],'.-',label=r"$\Lambda=100$")
axins.plot(strains4[:10],stresses4[:10],'.-',label=r"$\Lambda=1000$")

fig1.savefig("results/stress_strain.pdf")



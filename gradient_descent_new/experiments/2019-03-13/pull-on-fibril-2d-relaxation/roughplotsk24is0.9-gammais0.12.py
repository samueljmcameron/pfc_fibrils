import numpy as np
import matplotlib.pyplot as plt




data1 = np.loadtxt("data/_observables_scanforward_3.0000e+01_9.0000e-01_5.0000e+01_5.0000e+00_1.2000e-01.txt")
#data2 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+03_2.5000e+01_4.0000e-02.txt")
#data3 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+04_2.5000e+01_4.0000e-02.txt")


strains1 = data1[:,0]*100
energys1 = data1[:,1]
Rs1 = data1[:,2]
deltas1 = data1[:,4]
twists1 = data1[:,5]


#strains2 = data2[:,0]*100
#energys2 = data2[:,1]
#deltas2 = data2[:,4]
#twists2 = data2[:,5]


#strains3 = data3[:,0]*100
#energys3 = data3[:,1]
#deltas3 = data3[:,4]
#twists3 = data3[:,5]

fig,axarr = plt.subplots(4,sharex=True)
fig.set_size_inches(5,4*3.8)


axarr[0].plot(strains1,energys1,'o-',label=r"$\Lambda=50$")
#axarr[0].plot(strains2,energys2,'^-',label=r"$\Lambda=1000$")
#axarr[0].plot(strains3,energys3,'s-',label=r"$\Lambda=10000$")
axarr[0].set_ylabel(r"$E$")


axarr[1].plot(strains1,deltas1/deltas1[0],'o-',label=r"$\Lambda=50$")
#axarr[1].plot(strains2,deltas2/deltas2[0],'^-',label=r"$\Lambda=1000$")
#axarr[1].plot(strains3,deltas3/deltas3[0],'s-',label=r"$\Lambda=10000$")
axarr[1].set_ylabel(r"$\delta/\delta_0$")

axarr[2].plot(strains1,twists1,'o-',label=r"$\Lambda=50$")
#axarr[2].plot(strains2,twists2,'^-',label=r"$\Lambda=1000$")
#axarr[2].plot(strains3,twists3,'s-',label=r"$\Lambda=10000$")
axarr[2].set_ylabel(r"$\psi(R)$")

axarr[3].plot(strains1,Rs1,'o-',label=r"$\Lambda=50$")
#axarr[3].plot(strains2,twists2,'^-',label=r"$\Lambda=1000$")
#axarr[3].plot(strains3,twists3,'s-',label=r"$\Lambda=10000$")
axarr[3].set_ylabel(r"$R$")
axarr[3].set_xlabel("strain (%)")

axarr[1].legend(frameon=False,loc="lower left")


fig.subplots_adjust(left=0.2,right=0.95)
fig.savefig("results/roughplot.pdf")




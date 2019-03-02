import numpy as np
import matplotlib.pyplot as plt



data1 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+00_2.5000e+01_4.0000e-02.txt")
data2 = np.loadtxt("data/_observables_scanforward_3.0000e+01_1.0000e-01_1.0000e+01_2.5000e+01_4.0000e-02.txt")


strains1 = data1[:9,0]*100
energys1 = data1[:9,1]
deltas1 = data1[:9,4]
twists1 = data1[:9,5]


strains2 = data2[:9,0]*100
energys2 = data2[:9,1]
deltas2 = data2[:9,4]
twists2 = data2[:9,5]



fig,axarr = plt.subplots(3,sharex=True)
fig.set_size_inches(3.8,3*3.8)

axarr[0].plot(strains1,energys1,'o-',label=r"$\Lambda=1$")
axarr[0].plot(strains2,energys2,'^-',label=r"$\Lambda=10$")

axarr[0].set_ylabel(r"$E$")

axarr[1].plot(strains1,deltas1/deltas1[0],'o-',label=r"$\Lambda=1$")
axarr[1].plot(strains2,deltas2/deltas2[0],'^-',label=r"$\Lambda=10$")

axarr[1].set_ylabel(r"$\delta/\delta_0$")

axarr[2].plot(strains1,twists1,'o-',label=r"$\Lambda=1$")
axarr[2].plot(strains2,twists2,'^-',label=r"$\Lambda=10$")

axarr[2].set_ylabel(r"$\psi(R)$")
axarr[2].set_xlabel("strain (%)")

axarr[1].legend(frameon=False,loc="lower left")


fig.subplots_adjust(left=0.2,right=0.95)
fig.savefig("results/roughplot.pdf")




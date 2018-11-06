from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("QROMB_psivsr_3.0000e+01_8.0000e-01_1.0000e+00_1.0000e+00_1.0000e+00_8.7576e-02_6.3203e+00_7.9772e-01_5.0000e-02.txt")


r = data[:,0]
psi = data[:,1]
psiprime = data[:,2]
rf = data[:,3]

plt.plot(r,psi,'.')
plt.show()

plt.plot(r,psiprime,'.')
plt.show()

plt.plot(r,rf,'.')
plt.show()

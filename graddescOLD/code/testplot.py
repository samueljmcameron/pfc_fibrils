from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("../../../Data/Etest.txt")
data2 = np.loadtxt("../../../Data/EPsitest.txt")

plt.plot(data1[:,0],data1[:,1],'g.')

plt.figure()
plt.plot(data2[:,0],data2[:,1],'r.')

plt.figure()
plt.plot(data2[:,0],data2[:,2],'b.')

plt.figure()
plt.plot(data2[:,0],data2[:,2]**(3.0/5.0),'.')
plt.show()

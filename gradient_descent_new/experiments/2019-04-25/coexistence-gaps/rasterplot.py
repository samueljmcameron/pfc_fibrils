import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

def f(x,a,b,c):

    return a*x*x+b*x+c

if __name__ == "__main__":

    Lambda = 27
    omega = 10
    K33 = 30

    configure_fig_settings()

    data = np.loadtxt(f"data/_coexistence_{K33:1.4e}_{Lambda:1.4e}_"
                      +f"{omega:1.4e}.txt")


    gammas = data[:,0]
    k24s = data[:,1]


    popt = np.polyfit(k24s,gammas,2)

    k24s = np.linspace(0.5,1.0,num=51,endpoint=True)


    width = 3.37*2
    height = width

    fig = plt.figure()
    fig.set_size_inches(width,height)
    ax = fig.add_subplot(1,1,1)

    i = 0

    while (i < len(k24s)):

        k24 = k24s[i]
        gamma = f(k24,*popt)

        print(gamma,k24)

        gnews = np.linspace(gamma-0.002,gamma+0.002,num=7,endpoint=False)

        ax.plot(gnews,k24*np.ones(7),'.')


        i += 1

    ax.plot(f(k24s,*popt),k24s,'-')

    ax.set_xlabel(r"$\gamma$",fontsize=10)
    ax.set_ylabel(r"$k_{24}$",fontsize=10)


    fig.savefig("results/rasterplot-ex.pdf")

    plt.show()

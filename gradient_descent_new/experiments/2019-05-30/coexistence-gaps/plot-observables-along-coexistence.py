import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from readparams import ReadParams
import seaborn as sns

if __name__ == "__main__":

    Lambda = 27
    omega = 10
    K33 = 30

    width=3.37
    height=width/1.5

    configure_fig_settings()

    colors = sns.color_palette()

    colors = [colors[1],colors[2]] 

    scan = {}
    scan['\\Lambda'] = str(Lambda)
    scan['\\omega']= str(omega)
    scan['K_{33}'] = str(K33)

    loadsuf=["K_{33}","\\Lambda","\\omega"]
    savesuf=["K_{33}","\\Lambda","\\omega"]

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    suff = rp.write_suffix()

    fw_coexist = np.loadtxt(f"results/_fwd_coexist_{suff}.txt")
    bw_coexist = np.loadtxt(f"results/_bkwd_coexist_{suff}.txt")

    gammas = fw_coexist[:,0]

    observables_list = ["k24","E","R","eta","delta","psi(R)"]

    for i,observable in enumerate(observables_list):

        fig = plt.figure()
        fig.set_size_inches(width,height)

        ax = fig.add_subplot(1,1,1)
        print(observable)


        if i > 2:

            ylabel = f"\\{observable}"

        else:

            ylabel = observable

        if i == 0:

            ax.plot(gammas,fw_coexist[:,i+1],'.',color='k')
            ylabel="k_{24}"

        else:
            
            if observable == "psi(R)":

                linearlab = r"$\psi_{\ell}^*(R)$"
                frustlab = r"$\psi_{f}^*(R)$"

            else:

                linearlab = rf"${ylabel}_{{\ell}}^*$"
                frustlab = rf"${ylabel}_f^*$"

            ax.plot(gammas,fw_coexist[:,i+1],'^',color=colors[0],
                    label=linearlab,markersize=3)
            ax.plot(gammas,bw_coexist[:,i+1],'s',color=colors[1],
                    label=frustlab,markersize=3)

        if observable == "R":
            ax.set_yscale('log')

        ax.set_ylabel(rf"${ylabel}$",fontsize=10)
        ax.set_xlabel(rf"$\gamma$" " (along coexistence)",fontsize=10)
        ax.legend(frameon=False,fontsize=10)
        fig.subplots_adjust(left=0.2,bottom=0.2)
        
        fig.savefig(f"results/_{observable}_GAP_{suff}.pdf")




    

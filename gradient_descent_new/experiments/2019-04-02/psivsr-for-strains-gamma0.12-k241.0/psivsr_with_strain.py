import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from psidata import PsiData
from readparams import ReadParams

if __name__=="__main__":

    configure_fig_settings()

    width  = 3.37
    height = width


    k24 = '1.0'
    gamma = '0.12'
    omega = '10.0'
    Lambda = '100.0'

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24
    scan['\\omega']=omega
    scan['\\Lambda']=Lambda





    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = loadsuf


    #strains = np.linspace(0,0.034,num=35,endpoint=True)

    #strains = np.array([0,0.005,0.01,0.015],float)

    strains = np.array([0.0,0.01,0.02,0.03],float)

    colors = sns.color_palette("hls",len(strains))

    fig = plt.figure()
    fig.set_size_inches(width,height)
    
    ax = fig.add_subplot(1,1,1)



    obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)


    if obsfwd.E()[0] > 1e299:
        print("bad calculation at Lambda = ",Lambda)


    for i,u in enumerate(strains):

        if i == 0:

            strain = None

        else:
            
            strain = str(u)

            

        psi = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,strain=strain)

        ax.plot(psi.r(),psi.psi(),'-',color=colors[i],label=rf"$\epsilon={u:.3}$")

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\psi(r)$")
    ax.legend(frameon=False)
    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("psivsr-vsstrain",plot_format="pdf"))

    plt.show()









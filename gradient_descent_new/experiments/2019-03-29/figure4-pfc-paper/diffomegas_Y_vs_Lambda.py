import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from readparams import ReadParams


if __name__=="__main__":

    colors=sns.color_palette()
    
    configure_fig_settings()

    width  = 3.37
    height = width*1.25
    
    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i in range(2):
        
        ax[i] = fig.add_subplot(2,1,i+1)



    gammas = ['0.04','0.12']
    k24s = ['0.0','1.0']
    omegas = ['1.0','5.0','10.0','15.0','20.0','25.0','30.0']

    for ks,coord in enumerate(zip(gammas,k24s)):

        for ls,omega in enumerate(omegas):
        
            gamma = coord[0]
            k24 = coord[1]

            scan = {}
            scan['\\gamma_s']=gamma
            scan['k_{24}']=k24
            scan['\\omega']=omega

            Larray = np.linspace(0,1000,num=101,endpoint=True)
            Larray[0] = 1.0


            youngs = np.copy(Larray)*0

            Lambdas = list(map(str,Larray))

            loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
            savesuf = ["K_{33}","\\omega"]



            for js,Lambda in enumerate(Lambdas):

                scan['\\Lambda']=Lambda


                obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                        savesuf=savesuf)

                if len(obsfwd.data) == 6 or obsfwd.E()[0] > 1e299:
                    print("bad calculation at Lambda = ",Lambda)
                    Larray[js] = np.nan
                    continue
                strains = obsfwd.data[:,0]
                stresses = np.gradient(obsfwd.E(),strains)
                stresses[0] = 0.0
                youngs[js] = (stresses[1]-stresses[0])/(strains[1]-strains[0])


            ax[ks].plot(Larray,youngs,'o-',color=colors[ls],
                              label=rf"$\omega={omega}$")
            ax[ks].set_ylabel(r"$Y$")
            ax[ks].set_xscale('log')
            ax[ks].set_yscale('log')

            
    ax[1].set_xlabel(r"$\Lambda$")
    ax[1].set_xscale('log')
    ax[1].legend(loc="lower left",bbox_to_anchor=(0.2,0.7),
                        frameon=False)
    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("youngs-vs-Lambda-diffomega",plot_format="pdf"))

    plt.show()









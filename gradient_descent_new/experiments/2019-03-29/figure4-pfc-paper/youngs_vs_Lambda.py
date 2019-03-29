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

    configure_fig_settings()

    width  = 3.37
    height = width*1.25


    k24 = '1.0'
    gamma = '0.12'
    omega = sys.argv[1]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24
    scan['\\omega']=omega

    Larray = np.linspace(0,1000,num=101,endpoint=True)
    Larray[0] = 1.0


    youngs = np.copy(Larray)*0
    breaks = np.copy(Larray)*0

    colors = sns.color_palette("hls",len(Larray))

    Lambdas = list(map(str,Larray))

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]

    #observable_list = ["E","R","delta","surfacetwist","stress"]
    observable_list = ["youngs","breaks"]#,"delta","surfacetwist"]



    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i,observable in enumerate(observable_list):

        ax[observable] = fig.add_subplot(2,1,i+1)




    for js,Lambda in enumerate(Lambdas):

        scan['\\Lambda']=Lambda


        obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                savesuf=savesuf)
        #obsfwd.sort_observables()
        #obsfwd.remove_duplicates()

        if len(obsfwd.data) == 6 or obsfwd.E()[0] > 1e299:
            print("bad calculation at Lambda = ",Lambda)
            Larray[js] = np.nan
            continue
        strains = obsfwd.data[:,0]
        stresses = np.gradient(obsfwd.E(),strains)
        stresses[0] = 0.0
        youngs[js] = (stresses[1]-stresses[0])/(strains[1]-strains[0])

        breaks[js] = strains[np.argmax(stresses[stresses==stresses])]*100
        print(Lambda,breaks[js])        
        
    ax["youngs"].plot(Larray,youngs,'o-')
    ax["youngs"].set_ylabel(r"$Y$")
    ax["youngs"].set_xscale('log')
    ax["youngs"].set_yscale('log')
    ax["breaks"].plot(Larray,breaks,'o-')
    ax["breaks"].set_ylabel(r"$\epsilon_{\mathrm{max}}\times100\%$")
    ax["breaks"].set_xlabel(r"$\Lambda$")
    ax["breaks"].set_xscale('log')

    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("youngs-vs-Lambda",plot_format="pdf"))

    plt.show()









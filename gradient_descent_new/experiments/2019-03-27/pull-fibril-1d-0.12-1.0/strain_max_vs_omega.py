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
    height = width*1.5


    k24 = '1.0'
    gamma = '0.12'

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24

    omegaarray = np.array([1.0,5.0,10.0,15.0,20.0,25.0,30.0],float)


    colors = sns.color_palette("hls",len(omegaarray))

    omegas = list(map(str,omegaarray))

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\gamma_s"]

    observable_list = ["stress_max","strain_max"]


    ys_arrays = {}
    ys_arrays["strain_max"] = np.empty([len(omegas)],float)
    ys_arrays["stress_max"] = np.empty([len(omegas)],float)


    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i,observable in enumerate(observable_list):

        ax[observable] = fig.add_subplot(2,1,i+1)

    for kj,Lambda in enumerate(['1.0','10.0']):

        scan['\\Lambda']=Lambda

        for js,omega in enumerate(omegas):

            scan['\\omega']=omega
            print(omega)

            obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                    savesuf=savesuf)

            if not isinstance(obsfwd.E(),np.ndarray):

                ys_arrays["stress_max"][js]=0.0
                ys_arrays["strain_max"][js]=0.0
                continue

            if obsfwd.E()[0] > 1e299:
                print("bad calculation at Lambda = ",Lambda)
                continue

            strains = obsfwd.data[:,0]
            stresses = np.gradient(obsfwd.E(),strains)
            stresses[0] = 0.0


            max_index = np.argmax(stresses)


            ys_arrays["strain_max"][js] = strains[max_index]
            ys_arrays["stress_max"][js] = stresses[max_index]



        for i,observable in enumerate(observable_list):


            if observable == "stress_max":
                ylabel = r"$\sigma_{\mathrm{max}}$"
            elif observable == "strain_max":
                ylabel = r"$\epsilon(\sigma_{\mathrm{max}})$"


            ax[observable].plot(omegaarray,ys_arrays[observable],'.-',color=colors[kj],
                                    label=rf"$\Lambda={Lambda}$")
            ax[observable].set_ylabel(ylabel,fontsize = 10)
            ax[observable].set_xlabel(r"$\omega$",fontsize = 10)


            
    for observable in observable_list:

        ax[observable].set_xscale('log')
        if observable == "strain_max":
            ax[observable].legend(frameon=False)

    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("omega_dependent_failure",plot_format="pdf"))

    plt.show()






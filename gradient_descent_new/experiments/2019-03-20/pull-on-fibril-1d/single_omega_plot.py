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

    #configure_fig_settings()

    width  = 3.375
    height = width

    colors = sns.color_palette()


    k24 = '0.1'
    gamma = '0.04'
    omega = '10.0'

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24
    scan['\\omega']=omega

    Lambdas = ['1.0','10.0','100.0','1000.0']

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]

    observable_list = ["E","R","delta","surfacetwist","stress"]

    fig = {}
    ax = {}

    for i,observable in enumerate(observable_list):


        fig[observable] = plt.figure()
        ax[observable] = fig[observable].add_subplot(1,1,1)

        fig[observable].set_size_inches(width,height)


    for js,Lambda in enumerate(Lambdas):


        scan['\\Lambda']=Lambda

        obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                savesuf=savesuf)
        #obsfwd.sort_observables()
        #obsfwd.remove_duplicates()

        strains = obsfwd.data[:,0]
        stresses = np.gradient(obsfwd.E(),strains)
        stresses[0] = 0.0
        ysfwd = [obsfwd.E(),obsfwd.R(),obsfwd.delta(),obsfwd.surfacetwist(),stresses]



        for i,observable in enumerate(observable_list):


            if observable == 'surfacetwist':
                ylabel = r'$\psi(R)$'
            elif observable == "stress":
                ylabel = r"$\sigma$"
            elif observable == 'delta':
                ylabel = r'$\delta/\delta_0$'
                ysfwd[i] = ysfwd[i]/np.sqrt(2/3)
            elif len(observable) > 1:
                ylabel = fr'$\{observable}$'
            else:
                ylabel = fr'${observable}$'

            if observable == "stress":
                a = ysfwd[i][:]
                xs = strains[a>0]*100
                ys = a[a>0]
                ax[observable].plot(xs,ys,'.-',color=colors[js],
                                    label=rf"$\Lambda={Lambda}$")
            else:
                ax[observable].plot(strains*100,ysfwd[i][:],'.-',color=colors[js],
                                    label=rf"$\Lambda={Lambda}$")
            ax[observable].set_ylabel(ylabel,fontsize = 10)
            ax[observable].set_xlabel("strain (%)",fontsize = 10)

        #ax[observable].plot(Lambdas[:i_fwd+1],ysfwd[i][:i_fwd+1],'.',color=colors[i])
        #ax[observable].plot(Lambdas[i_fwd:i_fwd+2],ysfwd[i][i_fwd:i_fwd+2],'-',
        #                    color=colors[i])
        #ax[observable].plot(ms_bkwd,ysbkwd[i][1:],'-',color=colors[i])
        #ax[observable].plot(Lambdas[i_fwd+2:],ysfwd[i][i_fwd+2:],'.',color=colors[i])



        #ax[observable].tick_params("both",labelsize=18)
        

 
for observable in observable_list:

    if observable == "R":
        ax[observable].legend(frameon=False)
    if observable == "stress":
        ax[observable].set_yscale('log')
        ax[observable].set_xscale('log')    
    fig[observable].subplots_adjust(left=0.2,right=0.98,bottom=0.15,top=0.95)
    #fig[observable].savefig(obsfwd.observable_sname("gamma"+observable,plot_format="pdf"))



plt.show()






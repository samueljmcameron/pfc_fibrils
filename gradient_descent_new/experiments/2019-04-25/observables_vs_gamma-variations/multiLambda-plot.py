import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from psidata import PsiData
from readparams import ReadParams
from matplotlib.gridspec import GridSpec
from metastableregion import find_metastable

if __name__ == "__main__":

    k24 = sys.argv[1]
    omega = sys.argv[2]

    width = 3.37
    height = width

    configure_fig_settings()

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]
    savesuf = ["K_{33}","k_{24}","\\omega"]



    observable_list = ['E','R','eta','delta','surfacetwist']

    Lambdas = np.array([0.1,1.0,10.0,100.0],float)

    colors = sns.color_palette("hls",len(Lambdas))

    fig = {}

    ax = {}



    for observable in observable_list:

        fig[observable] = plt.figure()

        fig[observable].set_size_inches(width,height)

        ax[observable] = fig[observable].add_subplot(1,1,1)


    for j,Lambda in enumerate(Lambdas):

        scan = {}
        scan['\\Lambda']=Lambda
        scan['k_{24}']=k24
        scan['\\omega']=omega

        obsfwd = ObservableData(["\\gamma_s"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                savesuf=savesuf)

        gammas = obsfwd.data[:,0]

        ysfwd = [obsfwd.E(),obsfwd.R(),obsfwd.eta(),
                 obsfwd.delta(),obsfwd.surfacetwist()]

        """
        obsbkwd = ObservableData(["\\gamma_s"],scan_dir='scanbackward',scan=scan,loadsuf=loadsuf,
                                 savesuf=savesuf)

        ysbkwd = [obsbkwd.R(),obsbkwd.surfacetwist()]


        g_small,g_large,g_eq = find_metastable(obsfwd.R(),obsbkwd.R(),
                                               obsfwd.E(),obsbkwd.E())

        print(g_small,g_large,g_eq)
        """

        for i,observable in enumerate(observable_list):

            """
            ax[observable].plot(gammas[:g_eq],ysfwd[i][:g_eq],'-',color='k',
                                label=rf"$k_{{24}}={k24}$",
                                lw=3)
            ax[observable].plot(gammas[g_eq:g_large],ysfwd[i][g_eq:g_large],'-',color='k',
                                lw=1)
            ax[observable].plot(gammas[g_eq],ysfwd[i][g_eq],'^',color=colors[0],
                                markersize=3)
            ax[observable].plot(gammas[g_small:g_eq],ysbkwd[i][g_small:g_eq],'-',color='k',
                                lw=1)
            ax[observable].plot(gammas[g_eq:],ysbkwd[i][g_eq:],'-',color='k',
                                lw=3)
            ax[observable].plot(gammas[g_eq],ysbkwd[i][g_eq],'s',color=colors[1],
                                markersize=3)
            """

            ax[observable].plot(gammas,ysfwd[i][:],'.',color=colors[j],
                                label=rf"$\Lambda={Lambda}$")

            if observable == 'surfacetwist':
                ylabel = r'\psi(R)'
            elif len(observable)>1:
                ylabel = fr'\{observable}'
            else:
                ylabel = rf'{observable}'


            if observable == 'R':
                ax[observable].set_yscale('log')

            """
            ax[observable].annotate(linearlab,xy=(gammas[g_eq],ysfwd[i][g_eq]),
                                    xytext=(gammas[g_eq]+0.007,ysfwd[i][g_eq]*0.5),
                                    color=colors[0])#,
                                    #arrowprops=dict(facecolor='black', shrink=0.001))


            ax[observable].annotate(frustlab,xy=(gammas[g_eq],ysbkwd[i][g_eq]),
                                    xytext=(gammas[g_eq]-0.007,ysbkwd[i][g_eq]*1.3),
                                    color=colors[1])#,
                                    #arrowprops=dict(facecolor='black', shrink=0.001))
            """

            ax[observable].set_ylabel(rf"${ylabel}$",fontsize=10)
            ax[observable].set_xlabel(r"$\gamma$",fontsize = 10)
            ax[observable].legend(frameon=False,fontsize=10)

            ax[observable].set_xlim(0,0.3)


    for observable in observable_list:

        fig[observable].subplots_adjust(left=0.2,bottom=0.2)

        fig[observable].savefig(obsfwd.observable_sname("vary_Lambda"+observable,
                                                        plot_format="pdf"))



    plt.show()




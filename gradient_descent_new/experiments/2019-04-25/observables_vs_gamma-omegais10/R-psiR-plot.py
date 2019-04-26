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

omega = '10.0'

#height = 3.37/2
#width  = 4*height
#width = 3.37
#height = 5*width


width = 3.37
height = width/1.5

configure_fig_settings()

data2d = {}

colors = sns.color_palette()

colors = [colors[1],colors[2]]

loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]
savesuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]



observable_list = ['R','surfacetwist']

k24s = np.array([0.9],float)
Lambda = 27
omega = 10

fig = {}

ax = {}


for observable in observable_list:

    fig[observable] = plt.figure()
    
    fig[observable].set_size_inches(width,height)
    
    ax[observable] = fig[observable].add_subplot(1,1,1)


for j,k24 in enumerate(k24s):

    scan = {}
    scan['\\Lambda']=str(Lambda)
    scan['k_{24}']=str(k24)
    scan['\\omega']=str(omega)

    obsfwd = ObservableData(["\\gamma_s"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)

    gammas = obsfwd.data[:,0]

    ysfwd = [obsfwd.R(),obsfwd.surfacetwist()]


    obsbkwd = ObservableData(["\\gamma_s"],scan_dir='scanbackward',scan=scan,loadsuf=loadsuf,
                             savesuf=savesuf)

    ysbkwd = [obsbkwd.R(),obsbkwd.surfacetwist()]


    g_small,g_large,g_eq = find_metastable(obsfwd.R(),obsbkwd.R(),
                                           obsfwd.E(),obsbkwd.E())

    print(g_small,g_large,g_eq)

    for i,observable in enumerate(observable_list):


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


        if observable == 'surfacetwist':
            ylabel = r'\psi(R)'
            linearlab = r"$\psi_{\ell}^*(R)$"
            frustlab = r"$\psi_{f}^*(R)$"
        else:
            ylabel = fr'{observable}'

            linearlab = rf"${ylabel}_{{\ell}}^*$"
            frustlab = rf"${ylabel}_f^*$"


        if observable == 'R':
            ax[observable].set_yscale('log')

        ax[observable].annotate(linearlab,xy=(gammas[g_eq],ysfwd[i][g_eq]),
                                xytext=(gammas[g_eq]+0.007,ysfwd[i][g_eq]*0.5),
                                color=colors[0])#,
                                #arrowprops=dict(facecolor='black', shrink=0.001))


        ax[observable].annotate(frustlab,xy=(gammas[g_eq],ysbkwd[i][g_eq]),
                                xytext=(gammas[g_eq]-0.007,ysbkwd[i][g_eq]*1.3),
                                color=colors[1])#,
                                #arrowprops=dict(facecolor='black', shrink=0.001))


        ax[observable].set_ylabel(rf"${ylabel}$",fontsize=10)
        ax[observable].set_xlabel(r"$\gamma$",fontsize = 10)
        ax[observable].legend(frameon=False,fontsize=10)
        
        ax[observable].set_xlim(0,0.3)


for observable in observable_list:

    fig[observable].subplots_adjust(left=0.2,bottom=0.2)

    fig[observable].savefig(obsfwd.observable_sname(observable,plot_format="pdf"))



plt.show()




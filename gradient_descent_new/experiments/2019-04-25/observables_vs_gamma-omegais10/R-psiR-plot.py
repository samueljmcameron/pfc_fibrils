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
height = width

configure_fig_settings()

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]
savesuf = ["K_{33}","\\Lambda","\\omega"]



observable_list = ['R','surfacetwist']

k24s = np.array([1.0],float)
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


    g_small,g_large,g_eq = find_metastable(obsfwd,obsbkwd)

    print(g_small,g_large,g_eq)

    for i,observable in enumerate(observable_list):


        ax[observable].plot(gammas[:g_eq],ysfwd[i][:g_eq],'-',color=colors[j],
                            label=rf"$k_{{24}}={k24}$",
                            lw=3)
        ax[observable].plot(gammas[g_eq:g_large],ysfwd[i][g_eq:g_large],'-',color=colors[j],
                            label=rf"$k_{{24}}={k24}$",
                            lw=1)
        ax[observable].plot(gammas[g_small:g_eq],ysbkwd[i][g_small:g_eq],'-',color=colors[j],
                            label=rf"$k_{{24}}={k24}$",
                            lw=1)

        ax[observable].plot(gammas[g_eq:],ysbkwd[i][g_eq:],'-',color=colors[j],
                            label=rf"$k_{{24}}={k24}$",
                            lw=3)


        if observable == 'surfacetwist':
            ylabel = r'$\psi(R)$'
        else:
            ylabel = fr'${observable}$'

        if observable == 'R':
            ax[observable].set_yscale('log')

        ax[observable].set_ylabel(ylabel,fontsize=10)
        ax[observable].set_xlabel(r"$\gamma$",fontsize = 10)
    


for observable in observable_list:

    fig[observable].subplots_adjust(left=0.2,right=0.9,bottom=0.15,top=0.95)

    fig[observable].savefig(obsfwd.observable_sname(observable,plot_format="pdf"))



plt.show()




import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from readparams import ReadParams


gamma,k24,omega = '0.08','0.5','15'

width  = 8
height = width

configure_fig_settings()



fig = {}
ax = {}

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]
savesuf = loadsuf



Lambdas = np.linspace(0,1000,num=1001,endpoint=True)
i_fwd = 38
ms_fwd = np.array([38,39],float)
ms_bkwd = np.array([37,38,39,40],float)


scan = {}
scan['\\gamma_s']=gamma
scan['k_{24}']=k24    
scan['\\omega']=str(omega)

obsfwd = ObservableData(["\\Lambda"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                        savesuf=savesuf)
obsbkwd = ObservableData(["\\Lambda"],scan_dir='scanbackward',scan=scan,loadsuf=loadsuf,
                         savesuf=savesuf)
obsbkwd.sort_observables()

observable_list = ['E','R','eta','delta','surfacetwist']
ysfwd = [obsfwd.E(),obsfwd.R(),obsfwd.eta(),obsfwd.delta(),obsfwd.surfacetwist()]
ysbkwd = [obsbkwd.E(),obsbkwd.R(),obsbkwd.eta(),obsbkwd.delta(),obsbkwd.surfacetwist()]


for i,observable in enumerate(observable_list):


    fig[observable] = plt.figure()
    ax[observable] = fig[observable].add_subplot(1,1,1)

    fig[observable].set_size_inches(width,height)

    if observable == 'surfacetwist':
        ylabel = r'$\psi(R)$' + " (" + r"$^\circ$" + ")"
        ysfwd[i] *=180/np.pi
        ysbkwd[i] *= 180/np.pi
    elif observable == 'eta':
        ylabel = r'$2\pi/\eta$'
    elif observable == 'delta':
        ylabel = r'$\delta/\delta_0$'
    elif observable == 'R':
        ysfwd[i] /= 10 # divide by q (in units of inverse microns)
        ysfwd[i] *= 1000 # convert to nanometers
        ysbkwd[i] /= 10 # divide by q (in units of inverse microns)
        ysbkwd[i] *= 1000 # convert to nanometers
        ylabel = fr"${observable}$" + " (" + r"$\si{\nano\meter}$" + ")"

    elif len(observable) > 1:
        ylabel = fr'$\{observable}$'
    else:
        ylabel = fr'${observable}$'



    ax[observable].plot(Lambdas[:i_fwd+1],ysfwd[i][:i_fwd+1],'.',color=colors[i])
    ax[observable].plot(Lambdas[i_fwd:i_fwd+2],ysfwd[i][i_fwd:i_fwd+2],'-',
                        color=colors[i])
    ax[observable].plot(ms_bkwd,ysbkwd[i][1:],'-',color=colors[i])
    ax[observable].plot(Lambdas[i_fwd+2:],ysfwd[i][i_fwd+2:],'.',color=colors[i])
    ax[observable].set_xlabel(r"$\Lambda$",fontsize=24)
    ax[observable].set_ylabel(ylabel,fontsize=24)

    ax[observable].set_xscale('log')
    ax[observable].tick_params("both",labelsize=18)


 
for observable in observable_list:
    
    fig[observable].subplots_adjust(left=0.15,right=0.95,bottom=0.1)
    fig[observable].savefig(obsfwd.observable_sname(observable,plot_format="png"))

plt.show()




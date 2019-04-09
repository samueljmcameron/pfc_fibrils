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


omega = '10.0'

#height = 3.37/2
#width  = 4*height
#width = 3.37
#height = 5*width


width = 3.37#*2
height = 1.5*width

configure_fig_settings()



ax = {}

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]
savesuf = ["K_{33}","\\omega"]


#i_fwd = 38
#ms_fwd = np.array([38,39],float)
#ms_bkwd = np.array([37,38,39,40],float)


observable_list = ['E','eta','delta']


fig = plt.figure()

fig.set_size_inches(width,height)

for i,observable in enumerate(observable_list):

    ax[observable] = fig.add_subplot(3,1,i+1)




gammas = ['0.12']#['0.04','0.12']
k24s = ['1.0']#['0.0','1.0']

for js,pair in enumerate(zip(gammas,k24s)):

    js = js + 1

    if js == 0:
        Lambdalist=['184']
        mtypes = ["s"]
        ltypes = ["--"]
    elif js == 1:
        Lambdalist=['27']
        mtypes = ["X","D"]
        ltypes = ["-.",":"]
        ms_lower = 43
        ms_upper = 12

    i_Lambda = int(Lambdalist[0])

    gamma = pair[0]
    k24 = pair[1]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24    
    scan['\\omega']=str(omega)

    obsfwd = ObservableData(["\\Lambda"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)

    Lambdas = obsfwd.data[:,0]

    if js == 1:
        obsbkwd = ObservableData(["\\Lambda"],scan_dir='scanbackward',scan=scan,loadsuf=loadsuf,
                                 savesuf=savesuf)
        obsbkwd.sort_observables()

        ysbkwd = [obsbkwd.E(),obsbkwd.eta(),obsbkwd.delta()]



    ysfwd = [obsfwd.E(),obsfwd.eta(),obsfwd.delta()]


    for i,observable in enumerate(observable_list):


        if observable == 'surfacetwist':
            ylabel = r'$\psi(R)$'
        elif observable == 'eta':
            ylabel = r'$2\pi/\eta$'
            ysfwd[i] = 2*np.pi/ysfwd[i]
            if js == 1:
                ysbkwd[i] = 2*np.pi/ysbkwd[i]
        elif observable == 'delta':
            ylabel = r'$\delta/\delta_0$'
            ysfwd[i] = ysfwd[i]/np.sqrt(2/3)
            if js == 1:
                ysbkwd[i] = ysbkwd[i]/np.sqrt(2/3)
        elif len(observable) > 1:
            ylabel = fr'$\{observable}$'
        else:
            ylabel = fr'${observable}$'


        if js == 0:
            ax[observable].plot(Lambdas[1:],ysfwd[i][1:],'-',color=colors[js],
                                label=rf"$(\gamma,k_{{24}})=({gamma},{k24})$",
                                lw=3)
        else:
            ax[observable].plot(Lambdas[1:i_Lambda],ysfwd[i][1:i_Lambda],'-',color=colors[js],
                                label=rf"$(\gamma,k_{{24}})=({gamma},{k24})$",
                                lw=3)
            ax[observable].plot(Lambdas[i_Lambda:ms_lower],
                                ysfwd[i][i_Lambda:ms_lower],'-',color=colors[js],
                                lw=1)

            ax[observable].plot(Lambdas[ms_upper:i_Lambda],ysbkwd[i][ms_upper:i_Lambda],
                                '-',color=colors[js],lw=1)
            ax[observable].plot(Lambdas[i_Lambda:],ysbkwd[i][i_Lambda:],'-',color=colors[js],
                                lw=3)


        ax[observable].set_ylabel(ylabel,fontsize = 10)
        
    ax[observable].set_xlabel(r"$\Lambda$",fontsize = 10)
 
for observable in observable_list:

    if observable == "E":
        ax[observable].legend(frameon=False,loc="best",bbox_to_anchor=(0.2,-0.07,0.5,0.5))
    
    ax[observable].set_xscale('log')
fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
fig.savefig(obsfwd.observable_sname("E-eta-delta-twopair",plot_format="pdf"))



plt.show()




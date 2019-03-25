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

omega = '10.0'

#height = 3.37/2
#width  = 4*height
#width = 3.37
#height = 5*width


width = 3.37*2
height = width/2

configure_fig_settings()



ax = {}

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]
savesuf = ["K_{33}","\\omega"]
psi_loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]



#i_fwd = 38
#ms_fwd = np.array([38,39],float)
#ms_bkwd = np.array([37,38,39,40],float)


observable_list = ['R','surfacetwist']


fig = plt.figure()

fig.set_size_inches(width,height)

gs = GridSpec(2,5,figure=fig,hspace=0.5,wspace=1.5)

for i,observable in enumerate(observable_list):

    ax[observable] = fig.add_subplot(gs[i,:2])

ax[observable_list[0]].get_shared_x_axes().join(ax[observable_list[0]],
                                                ax[observable_list[1]])

ax['psi(r)'] = fig.add_subplot(gs[:,2:])


gammas = ['0.04','0.12']
k24s = ['0.0','1.0']

for js,pair in enumerate(zip(gammas,k24s)):

    gamma = pair[0]
    k24 = pair[1]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24    
    scan['\\omega']=str(omega)

    obsfwd = ObservableData(["\\Lambda"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)

    Lambdas = obsfwd.data[:,0]


    #obsbkwd = ObservableData(["\\Lambda"],scan_dir='scanbackward',scan=scan,loadsuf=loadsuf,
    #                         savesuf=savesuf)
    #obsbkwd.sort_observables()


    ysfwd = [obsfwd.R(),obsfwd.surfacetwist()]
    #ysbkwd = [obsbkwd.E(),obsbkwd.R(),obsbkwd.eta(),obsbkwd.delta(),obsbkwd.surfacetwist()]


    for i,observable in enumerate(observable_list):


        if observable == 'surfacetwist':
            ylabel = r'$\psi(R)$'
        elif observable == 'eta':
            ylabel = r'$2\pi/\eta$'
            ysfwd[i] = 2*np.pi/ysfwd[i]
        elif observable == 'delta':
            ylabel = r'$\delta/\delta_0$'
            ysfwd[i] = ysfwd[i]/np.sqrt(2/3)
        elif len(observable) > 1:
            ylabel = fr'$\{observable}$'
        else:
            ylabel = fr'${observable}$'


        ax[observable].plot(Lambdas[1:],ysfwd[i][1:],'.-',color=colors[js],
                            label=rf"$(\gamma,k_{{24}})=({gamma},{k24})$")

        #ax[observable].plot(Lambdas[:i_fwd+1],ysfwd[i][:i_fwd+1],'.',color=colors[i])
        #ax[observable].plot(Lambdas[i_fwd:i_fwd+2],ysfwd[i][i_fwd:i_fwd+2],'-',
        #                    color=colors[i])
        #ax[observable].plot(ms_bkwd,ysbkwd[i][1:],'-',color=colors[i])
        #ax[observable].plot(Lambdas[i_fwd+2:],ysfwd[i][i_fwd+2:],'.',color=colors[i])

        ax[observable].set_ylabel(ylabel,fontsize = 10)

        #ax[observable].tick_params("both",labelsize=18)
        
    ax[observable].set_xlabel(r"$\Lambda$",fontsize = 10)


    if js == 1:
        for Lambda in ['24.0']:# ['1.0','24.0','184.0','1000.0']:

            scan['\\Lambda'] = Lambda
            psistuff = PsiData(scan=scan,loadsuf=psi_loadsuf,savesuf=psi_loadsuf,
                               name=f"psivsr")
            ax["psi(r)"].plot(psistuff.r()/psistuff.r()[-1],
                              psistuff.psi(),label=rf"$\Lambda={Lambda}$",
                              color=colors[js])

            if Lambda == '24.0':

                psistuff = PsiData(scan=scan,loadsuf=psi_loadsuf,savesuf=psi_loadsuf,
                                   name=f"frustratedpsivsr")
                ax["psi(r)"].plot(psistuff.r()/psistuff.r()[-1],psistuff.psi(),
                                  label=rf"$\Lambda={Lambda}$",
                                  color=colors[js])

    

 
for observable in observable_list:

    if observable == "R":
        ax[observable].set_yscale('log')
        ax[observable].set_zorder(1)
        ax[observable].legend(frameon=True,loc="upper left",
                              bbox_to_anchor=(0.0,0.0))
        ax[observable].tick_params(labelbottom=False)
    
    ax[observable].set_xscale('log')


# now plot psi(r)

ax["psi(r)"].set_xlabel(r"$r/R$",fontsize=10)
ax["psi(r)"].set_ylabel(r"$\psi(r)$",fontsize=10)



 
fig.subplots_adjust(left=0.1,right=0.9,bottom=0.15,top=0.95)
fig.savefig(obsfwd.observable_sname("radial-twopair",plot_format="pdf"))



plt.show()




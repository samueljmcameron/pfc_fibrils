import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../modules_gammak24/')
from plotobservables import PlotObservables
from readparams import ReadParams


width  = 3.487
height = width

# see user_inputs.md for details on what typically goes in these inputs.
user_input = input("input string of a gamma,k24 pair, "
                   "using comma as delimiter: ")
gamma,k24 = user_input.split(',')


scan = {}
scan['\\gamma_s']=gamma
scan['k_{24}']=k24


observable_list = ['E','R','eta','delta','surfacetwist']

observable_levels = {}

Lambda_yticks = {}

omega_yticks = {}

observable_levels['E'] = np.array([-4,-3,-2,-1,0],float)

observable_levels['R'] = np.array([0.225,0.23,1.9,2.0],float)

observable_levels['eta'] = np.array([0.985,0.99,0.993],float)

observable_levels['delta'] = np.array([0.7,0.75,0.8,0.81],float)

observable_levels['surfacetwist'] = np.array([0.16,0.18,0.198],float)

Lambda_yticks['surfacetwist'] = [0.2,0.3]

omega_yticks['surfacetwist'] = [0.2,0.3]

Lambda_yticks['delta'] = [0.805,0.81,0.815]

omega_yticks['delta'] = [0,0.4,0.8]

Lambda_yticks['eta'] = [0.985,0.99,0.995]

omega_yticks['eta'] = [0.993,0.994]

Lambda_yticks['R'] = [1,2]

omega_yticks['R'] = [0.25,0.35]

Lambda_yticks['E'] = [-2.5,-2]

omega_yticks['E'] = [-4,-2,0]





configure_fig_settings()

fig = {}
ax_main = {}
ax_Lambda = {}
ax_omega = {}

data2d = {}

gs = gridspec.GridSpec(2,2,width_ratios=[3,1],height_ratios=[1,3])
gs.update(wspace=0.0,hspace=0.0)

for observable in observable_list:
    
    fig[observable] = plt.figure()

    fig[observable].set_size_inches(width,height)

    ax_main[observable] = plt.subplot(gs[1,0])
    ax_Lambda[observable] = plt.subplot(gs[0,0],
                                        sharex=ax_main[observable])
    ax_omega[observable] = plt.subplot(gs[1,1],
                                       sharey=ax_main[observable])


    data2d[observable]=np.empty([31,31],float)

colors = sns.color_palette()

savesuf = ["K_{33}","k_{24}","d_0","\\gamma_s"]
loadsuf = ["K_{33}","k_{24}","d_0","\\omega","\\gamma_s"]


num_Lambdas = 31
num_omegas = 31
Lambdas = np.linspace(0,30,num=num_Lambdas,endpoint=True)
omegas = np.linspace(0,30,num=num_omegas,endpoint=True)

for i,omega in enumerate(omegas):

    scan['\\omega']=str(omega)

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    obs = PlotObservables(["\\Lambda"],rp,scan_dir='scanforward')
    print(omega)
    for j,observable in enumerate(observable_list):
        
        index = obs.ylabelstr_to_column(observable)
        data2d[observable][i,:] = obs.data[:,index]
        if observable == 'eta':
            data2d[observable][i,:] = 2*np.pi/obs.data[:,index]
            data2d[observable][i,0] = np.nan



xlabel = r'$\Lambda$'
ylabel = r'$\omega$'



for observable in observable_list:

    if observable == 'surfacetwist':
        plabel = r'$\psi(R)$'
    elif observable == 'eta':
        plabel = r'$2\pi/\eta$'
    elif len(observable) > 1:
        plabel = fr'$\{observable}$'
    else:
        plabel = fr'${observable}$'



    z = data2d[observable]
    mask = np.zeros_like(z,dtype = bool)

    mask[:,0:2] = True
    mask[0:2,:] = True
    mask[:,24:25] = True
    mask[0:4,20:] = True
    mask[4:8,25:28] = True

    zmask = np.ma.array(z,mask=mask)

    cs = ax_main[observable].contourf(Lambdas,omegas,z)
    css = ax_main[observable].contour(Lambdas,omegas,zmask,
                                      observable_levels[observable],colors='r')
    ax_main[observable].clabel(css,fontsize=9,inline=1)

    cutout_omega_at_0 = data2d[observable][0,:]

    cutout_omega_at_15 = data2d[observable][15,:]

    cutout_omega_at_30 = data2d[observable][30,:]

    cutout_Lambda_at_0 = data2d[observable][:,0]

    cutout_Lambda_at_15 = data2d[observable][:,15]

    cutout_Lambda_at_30 = data2d[observable][:,30]


    ax_Lambda[observable].plot(Lambdas,cutout_omega_at_15,'.',
                               color='orange')
    plt.setp(ax_Lambda[observable].get_xticklabels(),visible=False)
    ax_Lambda[observable].set_yticks(Lambda_yticks[observable])
    ax_Lambda[observable].set_ylabel(plabel)
    ax_Lambda[observable].set_xlim(Lambdas[0],Lambdas[-1])
    ax_Lambda[observable].set_xticks([0,5,10,15,20,25,30])

    ax_main[observable].get_xticklabels()[3].set_color("magenta")
    
    ax_omega[observable].plot(cutout_Lambda_at_15,omegas,'.',
                              color='magenta')
    ax_omega[observable].set_xticks(omega_yticks[observable])
    ax_omega[observable].xaxis.tick_top()
    plt.setp(ax_omega[observable].get_yticklabels(),visible=False)
    ax_omega[observable].xaxis.set_label_position('top')
    ax_omega[observable].set_xlabel(plabel)
    ax_omega[observable].set_ylim(omegas[0],omegas[-1])
    
    ax_main[observable].get_yticklabels()[3].set_color("orange")

    ax_main[observable].set_xlabel(xlabel)
    ax_main[observable].set_ylabel(ylabel)
    fig[observable].subplots_adjust(left=0.18,right=0.95)
    fig[observable].savefig(obs.observable_sname(observable))

plt.show()




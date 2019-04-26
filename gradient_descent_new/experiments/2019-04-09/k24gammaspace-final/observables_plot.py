# for some reason, really tough to get this working well in manual contour labelling mode.


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from gridmbyn import GridmByn
from readparams import ReadParams
from matplotlib.colors import LinearSegmentedColormap

width  = 3.37
height = width

configure_fig_settings()

print('MANUAL CONTOURS ARE NOT WORKING WELL, SO BE SURE THAT THE FIGURES'
      'WITH THE CORRECT CONTOURS ARE SAVED BEFORE CONTINUING')


#observable_list = ['E','R','eta','delta','surfacetwist']
#observable_list = ['eta']

observable_levels = {}

gamma_yticks = {}

k24_yticks = {}

k24_ylims = {}

gamma_xlims = {}

vmins = {}

vmaxs = {}

vmins['E'] = -3

vmins['R'] = 0.01

vmins['surfacetwist'] = 0.05

vmins['eta'] = 0.8

vmins['delta'] = 0.95

vmaxs['E'] = -1

vmaxs['R'] = 100.0

vmaxs['surfacetwist'] = 0.60

vmaxs['eta'] = 1.0

vmaxs['delta'] = 1.0


k24_ylims['E'] = gamma_xlims['E'] = [-6,1]

k24_ylims['R'] = gamma_xlims['R'] = [0.02,5]

k24_ylims['surfacetwist'] = gamma_xlims['surfacetwist'] = [0,0.3]

k24_ylims['eta'] = gamma_xlims['eta'] = [0.97,1.0]

k24_ylims['delta'] = gamma_xlims['delta'] = [0.85,1.0]


observable_levels['E'] = np.array([-4,-2,-1.5,-1.22,-1.18,-1.17],float)

observable_levels['R'] = np.array([0.05,0.1,0.5,1.0,5.0,25.0],float)

observable_levels['eta'] = np.array([0.85,0.9,0.95,0.97,0.99,1.0],float)

observable_levels['delta'] = np.array([0.93,0.95,0.97,0.98,0.99],float)

observable_levels['surfacetwist'] = np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.5],float)

gamma_yticks['surfacetwist'] = [0.05,0.15,0.25]

k24_yticks['surfacetwist'] = [0.05,0.15,0.25]

gamma_yticks['delta'] = [0.9,0.95,1.0]

k24_yticks['delta'] = [0.9,0.95,1.0]

gamma_yticks['eta'] = [0.98,0.99,1.0]

k24_yticks['eta'] = [0.98,0.99,1.0]

gamma_yticks['R'] = [0.1,1.0]

k24_yticks['R'] = [0.1,1.0]

gamma_yticks['E'] = [-5,-1,3]

k24_yticks['E'] = [-5,-1,3]


fig = {}
ax = {}

data2d = {}

colors = sns.color_palette()

savesuf = ["K_{33}","\\Lambda","\\omega"]
loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]


num_gammas=101
gammas = np.linspace(0.01,0.4,num=101,endpoint=True)
num_k24s=101
k24s = np.linspace(-1,1,num=101,endpoint=True)



for observable in observable_list:
    
    fig[observable] = plt.figure()

    fig[observable].set_size_inches(width,height)

    ax[observable]= fig[observable].add_subplot(1,1,1)
    data2d[observable]=np.empty([num_k24s,num_gammas],float)




scan = {}
scan['\\Lambda']=sys.argv[1]
scan['\\omega']=sys.argv[2]



# load in 2d grid of data in data2d for each observable at the
# specified gamma,k24 pair.
for om,k24 in enumerate(k24s):
    
    scan['k_{24}']=str(k24)

    obs = ObservableData(["\\gamma_s"],scan=scan,loadsuf=loadsuf,savesuf=savesuf)
    print(k24)
    for observable in observable_list:

        index = obs.ylabelstr_to_column(observable)
        datalength = len(obs.data[:,index])
        if datalength < num_gammas-1:
            dumarray = np.ones(num_gammas-datalength)*np.nan
            dataarray = np.concatenate((obs.data[:,index],dumarray))
        else:
            dataarray = obs.data[:num_gammas,index]

        if observable == 'eta':
            data2d[observable][om,:] = 2*np.pi/dataarray
            if gammas[0] == 0.0:
                data2d[observable][om,0] = np.nan
        elif observable == 'delta':
            data2d[observable][om,:] = dataarray/np.sqrt(2/3)
        else:
            data2d[observable][om,:] = dataarray



xlabel = r'$\gamma$'
ylabel = r'$\k24$'

cmap = plt.get_cmap('Blues')
colors = cmap(np.linspace(0,0.75,cmap.N//2))
cmap2 = LinearSegmentedColormap.from_list('Upper Half',colors)


for observable in observable_list:

    if observable == 'surfacetwist':
        plabel = r'$\psi(R)$'
    elif observable == 'eta':
        plabel = r'$2\pi/\eta$'
    elif observable == 'delta':
        plabel = r'$\delta/\delta_0$'
    elif len(observable) > 1:
        plabel = fr'$\{observable}$'
    else:
        plabel = fr'${observable}$'

    z = data2d[observable]


    clev = np.linspace(vmins[observable],vmaxs[observable],num=101)


    cs = ax[observable].contourf(gammas,k24s,z,clev,cmap=cmap2,
                                 extend='both')

    css = ax[observable].contour(gammas,k24s,z,
                                 observable_levels[observable],
                                 colors='k',linestyles='dashed',linewidths=1)

    #ax[observable].clabel(css,manual=True,fontsize=10,inline=1,fmt='%1.2f')
    ax[observable].clabel(css,fontsize=10,inline=1,fmt='%1.2f')
    ax[observable].set_yticks([-1,-0.5,0,0.5,1])
    ax[observable].set_xlabel(r"$\gamma$",fontsize=10)
    ax[observable].set_ylabel(r"$k_{24}$",fontsize=10)
    

coexist = np.loadtxt(obs.observable_sname("coexistence",plot_format="txt"))

for observable in observable_list:
    
    ax[observable].plot(coexist[:,0],coexist[:,1],'k.')
    ax[observable].set_xlim(right=0.3)
    fig[observable].subplots_adjust(left=0.2,right=0.95,bottom=0.15)
    fig[observable].savefig(obs.observable_sname(observable))


plt.show()




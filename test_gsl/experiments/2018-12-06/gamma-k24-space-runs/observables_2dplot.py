import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../modules_gammak24/')
from plotobservables2d import PlotObservables2d


width  = 3.487
height = width

observable_list = ['energy','radius','eta','delta','surfacetwist']

configure_fig_settings()

fig = {}
ax = {}

for observable in observable_list:
    
    fig[observable],ax[observable] = plt.subplots()

    fig[observable].set_size_inches(width,height)


Lambda = float(sys.argv[1])
omega = float(sys.argv[2])
gamma_min,gamma_max = 0.01,0.2
k24_min,k24_max = 0.0,0.95


colors = sns.color_palette()

obs = PlotObservables2d(Lambda,omega,gamma_min,gamma_max,k24_min,k24_max)

xlabel=r'$\gamma$'
ylabel=r'$k_{24}$'

for j,observable in enumerate(observable_list):

    if observable == 'surfacetwist':
        obslabel = r'$\psi(R)$'
    elif observable == 'energy':
        obslabel = r'$E$'
    elif observable == 'radius':
        obslabel = r'$R$'
    elif len(observable) > 1:
        obslabel = fr'$\{observable}$'
    else:
        obslabel = fr'${observable}$'
    
    obs.load_data(observable)

    s = obs.plot_observable2d(ax[observable],observable,
                              label=obslabel)

    fig[observable].colorbar(s)
    fig[observable].suptitle(obslabel)
    ax[observable].set_xlabel(xlabel)
    ax[observable].set_ylabel(ylabel)
    fig[observable].tight_layout()
    fig[observable].savefig(obs.observable_sname(observable))



plt.show()

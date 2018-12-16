import numpy as np
import matplotlib.pyplot as plt
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

configure_fig_settings()

fig = {}
ax = {}

for observable in observable_list:
    
    fig[observable],ax[observable] = plt.subplots()

    fig[observable].set_size_inches(width,height)

colors = sns.color_palette()

savesuf = ["K_{33}","k_{24}","d_0","\\gamma_s"]
loadsuf = ["K_{33}","k_{24}","d_0","\\gamma_s"]



rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
obs = PlotObservables(["\\Lambda","\\omega"],rp)
print(obs.observables_fname())
for j,observable in enumerate(observable_list):
    obs.plot_observable(ax[observable],observable,color=colors[j],
                        label=fr'$\gamma_s,k_{{24}}={float(gamma):.2f},{float(k24):.1f}$')



xlabel = r'$3\Lambda=\omega$'



for observable in observable_list:

    if observable == 'surfacetwist':
        ylabel = r'$\psi(R)$'
    elif len(observable) > 1:
        ylabel = fr'$\{observable}$'
    else:
        ylabel = fr'${observable}$'

    ax[observable].set_xlabel(xlabel)
    ax[observable].set_ylabel(ylabel)
    ax[observable].legend(frameon=False)
    fig[observable].tight_layout()
    fig[observable].savefig(obs.observable_sname(observable))

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
user_input = input("input string of a Lambda,omega pair, "
                   "using comma as delimiter: ")
Lambda,omega = user_input.split(',')


scan = {}
scan['\\Lambda']=Lambda
scan['\\omega']=omega

observable_list = ['E','R','eta','delta','surfacetwist']

configure_fig_settings()

fig = {}
ax = {}

for observable in observable_list:
    
    fig[observable],ax[observable] = plt.subplots()

    fig[observable].set_size_inches(width,height)


k24s=np.linspace(0,.5,num=6,endpoint=True)

colors = sns.color_palette(n_colors=len(k24s))

savesuf = ["K_{33}","\\Lambda","d_0","\\omega"]
loadsuf = ["K_{33}","k_{24}","\\Lambda","d_0","\\omega"]



for i,k24 in enumerate(k24s):
    
    scan["k_{24}"] = float(k24)

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
    obs = PlotObservables("\\gamma_s",rp)
    print(obs.observables_fname())
    for j,observable in enumerate(observable_list):
        obs.plot_observable(ax[observable],observable,color=colors[i],
                            label=fr'$k_{{24}}={k24:.1f}$')



xlabel = r'$\gamma$'



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

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../local_packages/')
from plotobservables import PlotObservables



# see user_inputs.md for details on what typically goes in these inputs.
c_1 = input("please type the first gamma k24 coordinate, with a comma delimiter: ")

#c_2 = input("please type the second gamma k24 coordinate, with a comma delimiter: ")

#c_3 = input("please type the third gamma k24 coordinate, with a comma delimiter: ")

#c_list = [c_1,c_2,c_3]

c_list = [c_1]

coordinates = np.array([[float(s) for s in arg.split(',')] for arg in c_list],float)


width  = 3.487
height = width

observable_list = ['E','R','eta','delta','surfacetwist']

configure_fig_settings()

fig = {}
ax = {}

for observable in observable_list:
    
    fig[observable],ax[observable] = plt.subplots()

    fig[observable].set_size_inches(width,height)


for coordinate in coordinates:

    gamma,k24 = coordinate
    obs = PlotObservables(gamma,k24,scan_dir="scanbackward")
    obs.sort_observables()
    
    for observable in observable_list:
        obs.plot_observable_omega_eq_Lambda(ax[observable],observable,
                                            fr'$\gamma_s,k_{{24}}={gamma},{k24}$',
                                            'o-')

xlabel = r'$\omega=\Lambda$'

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



plt.show()

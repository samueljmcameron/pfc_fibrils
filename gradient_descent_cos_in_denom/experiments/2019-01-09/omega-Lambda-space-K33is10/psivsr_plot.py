
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../modules_gammak24/')
from plotpsi import PlotPsi
import seaborn as sns

colors = sns.color_palette()

configure_fig_settings()


# see user_inputs.md for details on what typically goes in these inputs.
user_input = input("input string of k24,Lambda,omega values, "
                   "using comma as delimiter: ")

gammas = np.array([0.025,0.05,0.075,0.1,0.125,0.15,0.175],float)


k24,Lambda,omega = np.array(user_input.split(','),float)


fig,ax = plt.subplots()
width  = 3.487
height = width
fig.set_size_inches(width,height)

for i,gamma in enumerate(gammas):
    if i == 0:
        label = fr'$\gamma=\num{{{gamma:.2e}}}$'
    elif i == len(gammas)-1:
        label = fr'$\gamma=\num{{{gamma:.2e}}}$'
    else:
        label = None
    pl = PlotPsi(gamma,k24,fig,ax)
    pl.plot_psivsr(omega,Lambda,label=label,color=colors[i])

ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$\psi(r)$')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(pl.psivsr_sname(omega,Lambda))

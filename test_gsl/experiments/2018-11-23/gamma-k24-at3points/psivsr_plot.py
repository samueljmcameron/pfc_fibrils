import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../local_packages/')
from plotpsi import PlotPsi

configure_fig_settings()


# see user_inputs.md for details on what typically goes in these inputs.
o_L_input = input("input string of omega Lambda values, "
                  "using comma as delimiter: ")

c_input = input("input string consisting of the gamma,k24 coordinate, "
                "using comma as delimiter: ")


omegaLambdas = np.array(o_L_input.split(','),float)
coordinate = np.array(c_input.split(','),float)


fig,ax = plt.subplots()
width  = 3.487
height = width
fig.set_size_inches(width,height)


pl = PlotPsi(*coordinate,fig,ax)

for ome_Lam in omegaLambdas:
    pl.plot_psivsr(ome_Lam,ome_Lam,
                   fr'$\omega=\Lambda=\num{{{ome_Lam:.1e}}}$')

ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$\psi(r)$')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(pl.psivsr_sname())

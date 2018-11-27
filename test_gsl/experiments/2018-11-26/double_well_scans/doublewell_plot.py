import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings

sys.path.append('../../local_packages/')
from plotdoublewell import PlotDoubleWell


if __name__=="__main__":

    configure_fig_settings()

    omegas_and_Lambdas = np.array([19.0,21.5,24.0],float)

    colors = sns.color_palette("muted",omegas_and_Lambdas.size)

    tmp = colors[1]
    colors[1] = colors[2]
    colors[2] = tmp
    
    fig,ax = plt.subplots()

    for i,o_L in enumerate(omegas_and_Lambdas):

        pdw = PlotDoubleWell(0.1,0.7,o_L,o_L,ax)

        pdw.plot_y_centred_dw(label=fr'$\omega=\Lambda=\num{{{o_L:.1f}}}$',
                              color=colors[i],markertype='.-')


    ax.legend(frameon=False)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$E$')
    ax.get_yaxis().set_ticks([])
    ax.get_xaxis().set_ticks([0,1])

    fig.savefig(pdw.sname())

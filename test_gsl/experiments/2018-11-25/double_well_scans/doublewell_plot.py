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

    omegas_and_Lambdas = np.linspace(19,24,num=6,endpoint=True)

    colors = sns.color_palette("muted",omegas_and_Lambdas.size)
    
    fig,ax = plt.subplots()

    for i,o_L in enumerate(omegas_and_Lambdas):

        pdw = PlotDoubleWell(0.1,0.7,o_L,o_L,ax)

        pdw.plot_y_centred_dw(label=fr'$\omega=\Lambda=\num{{{o_L:.0f}}}$',
                              color=colors[i])


    ax.legend(frameon=False)
    plt.show()

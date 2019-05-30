import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from plotdoublewell import PlotDoubleWell


if __name__=="__main__":

    pdw = PlotDoubleWell();

    ts = pdw.data[:,0]

    Es = pdw.data[:,1]

    configure_fig_settings()
    
    fig,ax = plt.subplots()


    width = 3.37
    height = 3.37

    fig.set_size_inches(width,height)


    ax.plot(ts,Es,'.')


    ax.legend(frameon=False)
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$E$')

    ax.get_xaxis().set_ticks([0,1])
    fig.subplots_adjust(left=0.3,bottom=0.2)

    fig.savefig(pdw.sname("Evst"))


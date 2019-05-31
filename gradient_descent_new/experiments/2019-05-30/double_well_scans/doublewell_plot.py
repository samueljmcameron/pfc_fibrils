import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from plotdoublewell import PlotDoubleWell


if __name__=="__main__":


    scan = {}

    scan['R0'] = str(4.29777048e-01)
    scan['R1'] = str(5.69248055e-01)
    scan['eta0'] = str(6.31908402e+00)
    scan['eta1'] = str(6.32575563e+00)
    scan['delta0'] = str(8.05969564e-01)
    scan['delta1'] = str(8.08987550e-01)
    scan['k_{24}'] = str(0.43)
    scan['\\gamma_s'] = str(8.39256250e-02)
    scan['\\Lambda'] = str(27.0)
    scan['\\omega'] = str(10.0)

    pdw = PlotDoubleWell(scan=scan);

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


import numpy as np
import matplotlib.pyplot as plt


class PlotPsi(object):

    def __init__(self,gamma,k24,fig,ax,K33=30.0,d0=1.0,lpath="data/",
                 spath="results/",plot_format="pdf"):

        self.gamma = gamma
        self.k24 = k24
        self.K33 = K33;
        self.d0 = d0;
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.fig = fig
        self.ax = ax

    def psivsr_fname(self,omega,Lambda):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{omega:.4e}"
                  f"_{self.d0:.4e}_{Lambda:.4e}_{self.gamma:.4e}")
    
        return f"{self.lpath}_psivsr{suffix}.txt"
    
    def psivsr_sname(self):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{self.d0:.4e}"
                  f"_{self.gamma:.4e}")
        return f"{self.spath}_psivsr{suffix}.{self.plot_format}"

    def plot_psivsr(self,omega,Lambda,label=None):

        data = np.loadtxt(self.psivsr_fname(omega,Lambda))
        rs = data[:,0]
        psis = data[:,1]
        self.ax.plot(rs,psis,'-',label=label)

        return

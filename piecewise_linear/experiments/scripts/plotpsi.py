# THIS IS DEPRECATED. EASIER TO JUST PLOT BY LOADING IN DATA USING psidata.py #


import numpy as np
import matplotlib.pyplot as plt


class PlotPsi(object):

    def __init__(self,gamma,k24,fig,ax,K33=30.0,lpath="data/",
                 spath="results/",plot_format="pdf"):

        self.gamma = gamma
        self.k24 = k24
        self.K33 = K33;
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.fig = fig
        self.ax = ax

    def psivsr_fname(self,omega,Lambda):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{omega:.4e}"
                  f"_{Lambda:.4e}_{self.gamma:.4e}")
    
        return f"{self.lpath}_psivsr{suffix}.txt"
    
    def psivsr_sname(self,omega,Lambda):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{omega:.4e}"
                  f"_{Lambda:.4e}")
        return f"{self.spath}_psivsr{suffix}.{self.plot_format}"

    def plot_psivsr(self,omega,Lambda,label=None,color=None):

        data = np.loadtxt(self.psivsr_fname(omega,Lambda))
        rs = data[:,0]
        psis = data[:,1]
        if color == None:
            self.ax.plot(rs,psis,'-',label=label)
        else:
            self.ax.plot(rs,psis,'-',label=label,color = color)

        return

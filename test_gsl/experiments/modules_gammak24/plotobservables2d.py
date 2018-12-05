import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class PlotObservables2d(object):

    def __init__(self,Lambda,omega,scan_dir="",data=None,K33=30.0,d0=1.0,lpath="data/",
                 spath="results/",plot_format="pdf",colors = sns.color_palette()):

        self.Lambda = Lambda
        self.omega = omega
        self.K33 = K33
        self.d0 = d0
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.scan_dir = f"_{scan_dir}"
        self.colors = colors
        self.data = data

    def load_data(self,varname):

        self.data = np.loadtxt(self.observables_fname(varname),dtype='float')

        return self.data

    def observables_fname(self,varname):

        suffix = (f"_{self.K33:.4e}_{self.Lambda:.4e}_{self.d0:.4e}"
                  f"_{self.omega:.4e}")
    
        return f"{self.lpath}_{varname}{suffix}.txt"
    
    def observable_sname(self,varname):

        suffix = (f"_{self.K33:.4e}_{self.Lambda:.4e}_{self.d0:.4e}"
                  f"_{self.omega:.4e}")

        return f"{self.spath}_{varname}{suffix}.{self.plot_format}"

    def plot_observable2d(self,ax,obs,label=None,switch_ysign=1):

        stuff = np.copy(self.data)
        """
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                if self.data[i,j] > 0:
                    stuff[i,j] = np.nan
                else:
                    stuff[i,j] = self.data[i,j]
                    """
        s=ax.imshow(stuff,interpolation='nearest',origin='lower')

        return s


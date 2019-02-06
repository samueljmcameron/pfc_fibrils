# THIS IS DEPRECATED. EASIER TO JUST PLOT BY LOADING IN DATA USING observabledata.py #


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class PlotObservables2d(object):

    def __init__(self,Lambda,omega,gamma_min,gamma_max,k24_min,k24_max,
                 scan_dir="",data=None,K33=30.0,lpath="data/",
                 spath="results/",plot_format="pdf",colors = sns.color_palette()):

        self.Lambda = Lambda
        self.omega = omega
        self.K33 = K33
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.scan_dir = f"_{scan_dir}"
        self.colors = colors
        self.data = data
        self.gamma_min = gamma_min
        self.gamma_max = gamma_max
        self.k24_min = k24_min
        self.k24_max = k24_max

    def load_data(self,varname):

        self.data = np.loadtxt(self.observables_fname(varname),dtype='float')

        return self.data

    def observables_fname(self,varname):

        suffix = (f"_{self.K33:.4e}_{self.Lambda:.4e}"
                  f"_{self.omega:.4e}")
    
        return f"{self.lpath}_{varname}{suffix}.txt"
    
    def observable_sname(self,varname):

        suffix = (f"_{self.K33:.4e}_{self.Lambda:.4e}"
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
        num_gamma = stuff.shape[0]
        num_k24 = stuff.shape[1]
        gammas = np.linspace(self.gamma_min,self.gamma_max,num=num_gamma,
                             endpoint=True)
        k24s = np.linspace(self.k24_min,self.k24_max,num=num_k24,
                           endpoint=True)

        aspect = (self.gamma_max-self.gamma_min)/(self.k24_max-self.k24_min)
        s=ax.imshow(stuff,interpolation='nearest',origin='lower',
                    extent = [self.gamma_min,self.gamma_max,self.k24_min,
                              self.k24_max],aspect=aspect)
        #gammas,k24s = np.meshgrid(gammas,k24s)
        CS = ax.contour(gammas,k24s,stuff,colors='k')
        ax.clabel(CS,fontsize=9,inline=1)
        return s


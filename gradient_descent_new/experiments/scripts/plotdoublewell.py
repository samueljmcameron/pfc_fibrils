import numpy as np
import matplotlib.pyplot as plt


class PlotDoubleWell(object):

    def __init__(self,gamma,k24,omega,Lambda,ax,K33=30.0,d0=1.0,lpath="data/",
                 spath="results/",plot_format="pdf"):

        self.gamma = gamma
        self.k24 = k24
        self.omega = omega
        self.Lambda = Lambda
        self.K33 = K33;
        self.d0 = d0;
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.ax = ax

    def fname(self):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{self.omega:.4e}"
                  f"_{self.Lambda:.4e}_{self.gamma:.4e}")
    
        return f"{self.lpath}_Evst{suffix}.txt"
    
    def sname(self):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}"
                  f"_{self.gamma:.4e}")
        return f"{self.spath}_Evst{suffix}.{self.plot_format}"

    def plot_dw(self,label=None,markertype='.',color=None):

        data = np.loadtxt(self.fname())
        ts = data[:,0]
        Es = data[:,1]
        if color == None:
            self.ax.plot(ts,Es,markertype,label=label)
        else:
            self.ax.plot(ts,Es,markertype,label=label,color=color)

        return


    def plot_y_centred_dw(self,label=None,markertype='.',color=None):

        data = np.loadtxt(self.fname())
        ts = data[:,0]
        Es = data[:,1]
        Es = Es-np.average(Es)
        if color == None:
            self.ax.plot(ts,Es,markertype,label=label)
        else:
            self.ax.plot(ts,Es,markertype,label=label,color=color)

        return

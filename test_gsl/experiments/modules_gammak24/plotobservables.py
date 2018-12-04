import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class PlotObservables(object):

    def __init__(self,gamma,k24,scan_dir="",data=None,K33=30.0,d0=1.0,lpath="data/",
                 spath="results/",plot_format="pdf",colors = sns.color_palette()):

        self.gamma = gamma
        self.k24 = k24
        self.K33 = K33;
        self.d0 = d0;
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.scan_dir = f"_{scan_dir}"
        self.colors = colors
        if (data==None):
            self.data = np.loadtxt(self.observables_fname())
        else:
            self.data = self.data

    def observables_fname(self):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{self.d0:.4e}"
                  f"_{self.gamma:.4e}")
    
        return f"{self.lpath}_observables{self.scan_dir}{suffix}.txt"
    
    def observable_sname(self,varname):

        suffix = f"_{self.K33:.4e}_{self.d0:.4e}"
        return f"{self.spath}_{varname}{suffix}.{self.plot_format}"

    def str_to_column(self,varname):
        
        if varname == 'omega':
            return 0
        elif varname == 'Lambda':
            return 1
        elif varname == 'E':
            return 2
        elif varname == 'R':
            return 3
        elif varname == 'eta':
            return 4
        elif varname == 'delta':
            return 5
        elif varname == 'surfacetwist':
            return 6
        else:
            raise ValueError(f"The label {varname} does not correspond to any "
                             f"columns in the file {self.observables_fname()}") 

    def plot_observable(self,ax,xvar_str,yvar_str,label=None,
                        markertype='.'):

        x = self.data[:,self.str_to_column(xvar_str)]
        y = self.data[:,self.str_to_column(yvar_str)]

        ax.plot(x,y,markertype,label=label)

        return

    def plot_observable_omega_eq_Lambda(self,ax,yvar_str,label=None,
                                        markertype='.',ms_markertype='-',
                                        color=None,switch_ysign=1,
                                        start_ms=None,end_ms=None):

        x = self.data[:,0]
        y = switch_ysign*self.data[:,self.str_to_column(yvar_str)]

        if start_ms==None:
            start_ms=len(x)
        if end_ms==None:
            end_ms=len(x)
        if color == None:
            color = self.colors[0]

        ax.plot(x[:start_ms],y[:start_ms],markertype,
                label=label,color = color)
        ax.plot(x[start_ms:end_ms],y[start_ms:end_ms],ms_markertype,
                color=color)
        ax.plot(x[end_ms:],y[end_ms:],markertype,color=color)


        return
    
    def sort_observables(self):

        self.data = self.data[self.data[:,0].argsort()]

        np.savetxt(self.observables_fname(),self.data,fmt='%13.6e',delimiter = '\t')

        return

# THIS IS DEPRECATED. EASIER TO JUST PLOT BY LOADING IN DATA USING observabledata.py #



import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from readparams import ReadParams

class PlotObservables(object):

    def __init__(self,xaxis,readparams,scan_dir="",data=None,
                 plot_format="pdf",colors = sns.color_palette()):

        self.xaxis=xaxis
        self.plot_format = plot_format
        self.scan_dir = scan_dir
        self.colors = colors
        self.readparams = readparams
        if (data==None):
            self.data = np.loadtxt(self.observables_fname())
        else:
            self.data = self.data

        return


    def observables_fname(self):
    
        suffix = self.readparams.write_suffix()

        return f"data/_observables_{self.scan_dir}_{suffix}.txt"
    
    def observable_sname(self,varname):

        suffix = self.readparams.write_suffix(suffix_type="save")

        return f"results/_{varname}_{suffix}.{self.plot_format}"

    def ylabelstr_to_column(self,varname,observables_num = 5):
        
        if isinstance(self.xaxis,list):
            vlength = len(self.xaxis)
        else:
            vlength = 1

        if varname == 'E':
            return vlength+0
        elif varname == 'R':
            return vlength+1
        elif varname == 'eta':
            return vlength+2
        elif varname == 'delta':
            return vlength+3
        elif varname == 'surfacetwist':
            return vlength+4
        else:
            raise ValueError(f"The label {varname} does not correspond to any "
                             f"columns in the file {self.observables_fname()}") 


    def plot_observable(self,ax,yvar_str,label=None,
                        markertype='.',ms_markertype='-',
                        color=None,switch_ysign=1,
                        start_ms=None,end_ms=None):

        x = self.data[:,0]
        y = self.data[:,self.ylabelstr_to_column(yvar_str)]

        y *= switch_ysign


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

# hello
# new version of plotobservables.py

# this reads in data, which can be accessed easily by calling e.g.
# ObservableData.E() for energy array, ObservableData.R() for radius
# array, etc.

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from readparams import ReadParams

class ObservableData(ReadParams):

    def __init__(self,xaxis=None,scan_dir="",datfile="data/input.dat",scan={},
                 loadsuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"]):

        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.xaxis = xaxis
        self.scan_dir = scan_dir
        self.data = np.loadtxt(self.observables_fname())

        return


    def observables_fname(self):
    
        suffix = self.write_suffix()

        if self.scan_dir != "":
            fname =f"data/_observables_{self.scan_dir}_{suffix}.txt"
        else:
            fname =f"data/_observables_{suffix}.txt"
        return fname
    
    def observable_sname(self,varname):

        suffix = self.write_suffix(suffix_type="save")

        return f"results/_{varname}_{suffix}.{self.plot_format}"

    def ylabelstr_to_column(self,varname,observables_num = 5):
        
        if isinstance(self.xaxis,list):
            vlength = len(self.xaxis)
        elif self.xaxis == None:
            vlength=0
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

    def E(self,observables_num=5):


        i = self.ylabelstr_to_column('E',observables_num=observables_num)

        dim = len(self.data.shape)
        
        if dim == 1:
            result = self.data[i]
        else:
            result = self.data[:,i]

        return result

    def R(self,observables_num=5):
        
        i = self.ylabelstr_to_column('R',observables_num=observables_num)

        dim = len(self.data.shape)
        
        if dim == 1:
            result = self.data[i]
        else:
            result = self.data[:,i]

        return result

    def eta(self,observables_num=5):
        
        i = self.ylabelstr_to_column('eta',observables_num=observables_num)

        dim = len(self.data.shape)
        
        if dim == 1:
            result = self.data[i]
        else:
            result = self.data[:,i]

        return result

    def delta(self,observables_num=5):
        
        i = self.ylabelstr_to_column('delta',observables_num=observables_num)

        dim = len(self.data.shape)
        
        if dim == 1:
            result = self.data[i]
        else:
            result = self.data[:,i]

        return result

    def surfacetwist(self,observables_num=5):

        i = self.ylabelstr_to_column('surfacetwist',observables_num=observables_num)

        dim = len(self.data.shape)
        
        if dim == 1:
            result = self.data[i]
        else:
            result = self.data[:,i]

        return result
    
    def sort_observables(self):

        self.data = self.data[self.data[:,0].argsort()]

        np.savetxt(self.observables_fname(),self.data,fmt='%13.6e',delimiter = '\t')

        return

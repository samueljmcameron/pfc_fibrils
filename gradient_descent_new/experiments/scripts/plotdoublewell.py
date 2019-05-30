import numpy as np
import matplotlib.pyplot as plt
from readparams import ReadParams
import os


class PlotDoubleWell(ReadParams):

    def __init__(self,xaxis=None,scan_dir="",datfile="data/input.dat",scan={},
                 loadsuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 name= "Evst"):

        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.xaxis = xaxis
        self.scan_dir = scan_dir
        self.name = name
        if os.path.isfile(self.fname()):
            self.data = np.loadtxt(self.fname())
            self.file_exists = True
        else:
            print("could not find a file by the name of ",
                  self.fname())
            self.file_exists = False

        return

    def fname(self):
    
        suffix = self.write_suffix()

        if self.scan_dir != "":
            fname =f"data/_{self.name}_{self.scan_dir}_{suffix}.txt"
        else:
            fname =f"data/_{self.name}_{suffix}.txt"
        return fname
    
    def sname(self,varname,plot_format="pdf"):

        suffix = self.write_suffix(suffix_type="save")

        return f"results/_{varname}_{suffix}.{plot_format}"


# new version of plotpsi.py

# this reads in data, which can be accessed easily by calling e.g.
# PsiData.r() for radial distance array, PsiData.psi() for psi(r)
# array, etc.

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from readparams import ReadParams


class PsiData(ReadParams):

    def __init__(self,scan_dir="",datfile="data/input.dat",scan={},
                 loadsuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],sfile_format=".pdf",
                 name="psivsr",strain=None):

        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.sfile_format = sfile_format
        self.scan_dir = scan_dir
        self.name = name
        self.strain = strain
        self.data = np.loadtxt(self.psivsr_fname())

        return

    def r(self):
        return self.data[:,0]

    def psi(self):
        return self.data[:,1]

    def psiprime(self):
        return self.data[:,2]

    def rf_fibril(self):
        return self.data[:,3]

    def psivsr_fname(self):

        suffix = self.write_suffix()

        if self.strain != None:

            strain = float(self.strain)

            suffix = suffix + f"_{strain:.4e}"
    
        return f"data/_{self.name}_{suffix}.txt"
    
    def psivsr_sname(self):

        suffix = self.write_suffix(suffix_type="save")
        
        if self.strain != None:

            suffix = suffix + f"_{self.strain:.4e}"

        if self.scan_dir != "":
            sname = f"results/_{self.name}_{self.scan_dir}_{suffix}.{self.sfile_format}"
        else:
            sname = f"results/_{self.name}_{suffix}.{self.sfile_format}"

        return sname




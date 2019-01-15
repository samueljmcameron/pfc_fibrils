import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from psidata import PsiData
from observabledata import ObservableData

class PlotFibrilStress(object):

    def __init__(self,psidata,observabledata,scan_dir="",
                 sfile_format=".pdf"):
        
        self.psidata = psidata
        self.observabledata = observabledata
        self.scan_dir=scan_dir
        self.sfile_format=sfile_format
        return


    def stresses(self):

        return 2*np.pi/np.cos(self.psidata.psi())-self.observabledata.eta()

    def plot(self,ax):

        azm = np.linspace(0,2*np.pi)

        rs,thetas=np.meshgrid(self.psidata.r(),azm)

        stresses= np.tile(self.stresses(),(rs.shape[0],1))

        im = ax.contourf(thetas,rs,stresses,cmap='bwr')

        return im

    def stress_sname(self):

        suffix = self.psidata.write_suffix(suffix_type="save")

        if self.scan_dir != "":
            sname = f"results/_stress_{self.scan_dir}_{suffix}.{self.sfile_format}"
        else:
            sname = f"results/_stress_{suffix}.{self.sfile_format}"

        return sname

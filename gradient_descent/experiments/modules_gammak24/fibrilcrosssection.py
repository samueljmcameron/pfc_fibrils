import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from psidata import PsiData
from observabledata import ObservableData

class FibrilStrain(object):

    def __init__(self,psidata,observabledata,scan_dir="",
                 sfile_format=".pdf"):
        
        self.psidata = psidata
        self.observabledata = observabledata
        self.scan_dir=scan_dir
        self.sfile_format=sfile_format
        return


    def mesh_polar(self,num_azm=50):

        azm = np.linspace(0,2*np.pi,num=num_azm)

        rs,thetas = np.meshgrid(self.psidata.r(),azm)

        return rs,thetas

    def strain_1d(self,denom='d(r)'):

        preferred_dband = np.cos(self.psidata.psi())
        true_dband = 2*np.pi/self.observabledata.eta()


        if denom=='d(r)':
            dn = preferred_dband
        elif denom=='d':
            dn = true_dband
        else:
            raise ValueError(f'setting denominator of the strain equation to '
                             '"{denom}" is not valid, it must either be "d(r)" '
                             '(the default value) or "d".')

        return (true_dband-preferred_dband)/dn

    def strain_polar(self,r_mesh,denom='d(r)'):

        return np.tile(self.strain_1d(denom=denom),
                       (r_mesh.shape[0],1))


    def strain_sname(self,descriptor="polar"):

        suffix = self.psidata.write_suffix(suffix_type="save")

        if self.scan_dir != "":
            sname = f"results/_{descriptor}_strain_{self.scan_dir}_{suffix}.{self.sfile_format}"
        else:
            sname = f"results/_{descriptor}_strain_{suffix}.{self.sfile_format}"

        return sname

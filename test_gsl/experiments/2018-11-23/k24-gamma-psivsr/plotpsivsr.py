import numpy as np
import matplotlib.pyplot as plt


args = ['fname']

args.append('0.15,0.5')
args.append('0.1,0.7')
args.append('0.02,0.9')

arg2 = ['0,0.1,1.0,10.0,100.0,1000.0']

a = np.array([arg.split(',') for arg in args[1:]],float)

omega_and_Lambda = np.array([arg.split(',') for arg in arg2][0],float)

print(omega_and_Lambda)


class PlotPsi(object):

    def __init__(self,gamma,k24,fig,ax,K33=30.0,d0=1.0,lpath="data/",
                 spath="results/",plot_format="pdf"):

        self.gamma = gamma
        self.k24 = k24
        self.K33 = K33;
        self.d0 = d0;
        self.spath = spath
        self.lpath = lpath
        self.plot_format = plot_format
        self.fig = fig
        self.ax = ax

    def psivsr_fname(self,omega,Lambda):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{omega:.4e}"
                  f"_{self.d0:.4e}_{Lambda:.4e}_{self.gamma:.4e}")
    
        return f"{self.lpath}_psivsr{suffix}.txt"
    
    def psivsr_sname(self):

        suffix = (f"_{self.K33:.4e}_{self.k24:.4e}_{self.d0:.4e}"
                  f"_{self.gamma:.4e}")
        return f"{self.spath}_psivsr{suffix}.{self.plot_format}"

    def plot_psivsr(self,omega,Lambda):

        data = np.loadtxt(self.psivsr_fname(omega,Lambda))
        rs = data[:,0]
        psis = data[:,1]
        self.ax.plot(rs,psis,'-')

        return


fig,ax = plt.subplots()
pl = PlotPsi(0.15,0.5,fig,ax)


for ome_Lam in omega_and_Lambda:
    pl.plot_psivsr(ome_Lam,ome_Lam)

plt.show()

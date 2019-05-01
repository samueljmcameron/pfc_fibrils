import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from psidata import PsiData
from observabledata import ObservableData
from fibrilstrain import FibrilStrain
from midpointnormalize import MidpointNormalize
import seaborn as sns
from matplotlib.gridspec import GridSpec

colors = sns.color_palette()

configure_fig_settings()



if len(sys.argv) < 5:

    user_input = input("input string of gamma,k24,Lambda,omega values, "
                       "using comma as delimiter: ")
    gamma,k24,Lambda,omega = np.array(user_input.split(','),float)
    denom = input('Now input either "d(r)" or "d" to set the denominator of the strain:\n')

else:

    gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]


width  = 3.37
height = width*2/3.0

fig = plt.figure()
fig.set_size_inches(width,height)
gs = GridSpec(1,3,figure=fig,wspace=0.0)
axarr = []

axarr.append(fig.add_subplot(gs[0,0]))
axarr.append(fig.add_subplot(gs[0,1:],projection='polar'))

scan = {}
scan['\\gamma_s'] = gamma
scan['k_{24}'] = k24
scan['\\Lambda'] = Lambda
scan['\\omega'] = omega


loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                   name=f"psivsr")


rs = psistuff.r()
psis = psistuff.psi()

observablestuff = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                                 name=f"observables")

R = observablestuff.R()
eta = observablestuff.eta()

print("radius is ",R)

dband = 2*np.pi/eta

strain_r = (dband-np.cos(psis))/np.cos(psis)

axarr[0].plot(rs,strain_r*100,'-',
              color=colors[0])

print("tension at fibril centre is ",strain_r[0]*100, "%.")
print("tension at the fibril surface is ",strain_r[-1]*100,"%.")


axarr[0].set_xlabel(r'$r$',fontsize=10)
axarr[0].set_ylabel('molecular strain (\%)',fontsize=10)
axarr[0].set_xlim(left=0)
#axarr[0].legend(frameon=False)

fibrilstrain = FibrilStrain(psistuff,observablestuff,sfile_format='pdf')

rs,thetas = fibrilstrain.mesh_polar(grid_skip=4)

strains = fibrilstrain.strain_polar(rs,grid_skip=4)

norm = MidpointNormalize(midpoint=0)

im = axarr[1].contourf(thetas,rs,strains,100,norm=norm,
                       cmap='bwr')

#clb = fig.colorbar(im,ax=axarr[i+1])

#clb.ax.set_title(rf'$\frac{{d-d(r)}}{{{denom}}}$')

axarr[1].set_xticks([])
axarr[1].set_yticks([])


#axarr[1].annotate(rf'$R={R:1.3f}$',xy=(5*np.pi/4,R-0.01*R),
#                  xytext=(5*np.pi/4,R+R-0.01*R),fontsize=20,
#                  color=colors[0])


fig.subplots_adjust(left=0.2,bottom=0.2,right=0.95)

fig.savefig(fibrilstrain.strain_sname(descriptor=f'polar'))


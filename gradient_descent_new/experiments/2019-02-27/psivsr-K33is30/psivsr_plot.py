
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


colors = sns.color_palette()

configure_fig_settings()



if len(sys.argv) < 6:

    user_input = input("input string of gamma,k24,Lambda,omega values, "
                       "using comma as delimiter: ")
    gamma,k24,Lambda,omega = np.array(user_input.split(','),float)
    denom = input('Now input either "d(r)" or "d" to set the denominator of the strain')
else:

    gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
    denom = sys.argv[5]


scan = {}
scan['\\gamma_s'] = gamma
scan['k_{24}'] = k24
scan['\\Lambda'] = Lambda
scan['\\omega'] = omega


loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)

rs = psistuff.r()/8.0*1000
psis = psistuff.psi()

observablestuff = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)

print(observablestuff.surfacetwist())

R = observablestuff.R()/8.0*1000


fig = plt.figure()
width  = 3.487
height = width
fig.set_size_inches(2*width,height)

ax1 = fig.add_subplot(1,2,1)



ax1.plot(rs,psis,'.',label=rf'$\Lambda={Lambda}$')

ax1.set_xlabel(r'$r$' ' ' r'$(\si{\nano\meter})$',fontsize=24)
ax1.set_ylabel(r'$\psi(r)$',fontsize=24)
ax1.legend(frameon=False)



ax2 = fig.add_subplot(1,2,2,projection='polar')

fibrilstrain = FibrilStrain(psistuff,observablestuff,sfile_format='.pdf')

rs,thetas = fibrilstrain.mesh_polar()

rs = rs/8.0*1000

strains = fibrilstrain.strain_polar(rs,denom=denom)

norm = MidpointNormalize(midpoint=0)

im = ax2.contourf(thetas,rs,strains,100,norm=norm,
                  cmap='bwr')

#clb = fig.colorbar(im,ax=ax2)

#clb.ax.set_title(rf'$\frac{{d-d(r)}}{{{denom}}}$')

ax2.set_xticks([])
ax2.set_yticks([])


ax2.annotate(rf'$R={R:1.0f}\si{{\nano\meter}}$',xy=(5*np.pi/4,R-0.01*R),xytext=(5*np.pi/4,R+2/3*R-0.01*R),fontsize=20)


fig.subplots_adjust(left=0.1,bottom=0.2)

fig.savefig(fibrilstrain.strain_sname(descriptor=f'polar_{denom}'))


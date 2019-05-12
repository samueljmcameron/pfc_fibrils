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

gamma = sys.argv[1]

k24 = sys.argv[2]
Lambda = sys.argv[3]
omega = sys.argv[4]

fig = plt.figure()
width  = 3.37
height = width
fig.set_size_inches(width,height)

ax1 = fig.add_subplot(1,1,1)

loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
savesuf=["K_{33}","k_{24}","\\Lambda","\\omega"]


q = 110/1000 # nm^{-1}



scan = {}
scan['\\gamma_s'] = gamma
scan['k_{24}'] = k24
scan['\\Lambda'] = Lambda
scan['\\omega'] = omega





psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,name=f"psivsr",
                   sfile_format="pdf")

rs = psistuff.r()
psis = psistuff.psi()

ax1.plot(rs/q,psis,'k-',label=rf'$\gamma={gamma}$',lw=2)

ax1.set_xlabel(r'$r$' + ' (' + r'$\si{\nano\meter}$' + ')',fontsize=10)
ax1.set_ylabel(r'$\psi(r)$' + ' (' + r'$\si{\radian}$' + ')',fontsize=10)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.legend(frameon=False,fontsize=10)


fig.subplots_adjust(left=0.2,bottom=0.2)

fig.savefig(psistuff.psivsr_sname())


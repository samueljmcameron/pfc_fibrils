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



if len(sys.argv) < 5:

    user_input = input("input string of gamma,k24,Lambda,omega values, "
                       "using comma as delimiter: ")
    gamma,k24,Lambda,omega = np.array(user_input.split(','),float)
else:

    gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

type = input("specify type (either linear, constant, "
             "or no type in which case just press enter):\n")


scan = {}
scan['\\gamma_s'] = gamma
scan['k_{24}'] = k24
scan['\\Lambda'] = Lambda
scan['\\omega'] = omega


R_units = 1000.0/10.0  # units of nano meters, with q = 10 (um)^{-1}

loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,name=f"psivsr{type}",
                   sfile_format=".png")

rs = psistuff.r()*R_units
psis = psistuff.psi()*180/np.pi



fig = plt.figure()
width  = 5
height = width
fig.set_size_inches(width,height)

ax1 = fig.add_subplot(1,1,1)



ax1.plot(rs,psis,'.',label=rf'$\Lambda={Lambda}$')

ax1.set_xlabel(r'$r$' ' ' r'$(\si{\nano\meter})$',fontsize=24)
ax1.set_ylabel(r'$\psi(r)$' + " (" + r"$^\circ$" + ")",fontsize=24)
ax1.legend(frameon=False,fontsize=24)
ax1.tick_params("both",labelsize=18)


fig.subplots_adjust(left=0.15,bottom=0.2)

fig.savefig(psistuff.psivsr_sname())


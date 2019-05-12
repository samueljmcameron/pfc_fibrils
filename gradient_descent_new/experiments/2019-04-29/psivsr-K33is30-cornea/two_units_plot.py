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

gamma = 0.15

k24s = [0.8,-0.3]
Lambdas = [0.9,1.5]
omega = 20

fig = plt.figure()
width  = 3.37
height = width
fig.set_size_inches(width,height)

ax1 = fig.add_subplot(1,1,1)

loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
savesuf = ["K_{33}","\\omega","\\gamma_s"]

markertypes=['--','-.']


qs = [24.0/1000,150.0/1000] # nm ^{-1}

types = ['linear','frustrated']

for i,type in enumerate(types):

    k24 = k24s[i]
    Lambda = Lambdas[i]
    q = qs[i]

    scan = {}
    scan['\\gamma_s'] = gamma
    scan['k_{24}'] = k24
    scan['\\Lambda'] = Lambda
    scan['\\omega'] = omega



    print(f"k24 = {k24}, Lambda = {Lambda}, q = {q}")


    psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,name=f"psivsr",
                       sfile_format="pdf")

    rs = psistuff.r()/q
    psis = psistuff.psi()

    ax1.plot(rs,psis,markertypes[i],label=f'{type} cornea',lw=2)

    print("Radius = ",rs[-1], " nm.")

    print("surface twist = ",psis[-1]," rad.")

ax1.set_xlabel(r'$r$' + ' (' + r'$\si{\nano\meter}$' + ')',fontsize=10)
ax1.set_ylabel(r'$\psi(r)$' + ' (' + r'$\si{\radian}$' + ')',fontsize=10)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.legend(frameon=False,fontsize=10)


fig.subplots_adjust(left=0.2,bottom=0.2)

psistuff.name = f"twocorneas"

fig.savefig(psistuff.psivsr_sname())



import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
import seaborn as sns


colors = sns.color_palette()

configure_fig_settings()



fig = plt.figure()
width  = 4
height = width
fig.set_size_inches(width,height)

ax1 = fig.add_subplot(1,1,1)


rs = np.linspace(0,0.1,num=101,endpoint=True)
psis = 1.2*rs



ax1.plot(rs,psis,'-',color=colors[0],lw=4)

ax1.set_xlabel(r'$r$',fontsize=24)
ax1.set_ylabel(r'$\psi(r)$',fontsize=24)
ax1.set_xticks([0,rs[-1]])
ax1.set_xticklabels([0,r"$R$"])
ax1.set_yticks([0,psis[-1]])
ax1.set_yticklabels([0,r"$\psi(R)$"])
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)


fig.subplots_adjust(left=0.2,bottom=0.15)

fig.savefig('results/constant_pitch.pdf')


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from psidata import PsiData
from readparams import ReadParams
from matplotlib.gridspec import GridSpec

omega = '10.0'

width = 3.37
height = width*2/3

configure_fig_settings()



ax = {}

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]
savesuf = ["K_{33}","\\omega"]
psi_loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]



observable_list = ['surfacetwist','eta']


fig = plt.figure()

fig.set_size_inches(width,height)

ax = fig.add_subplot(1,1,1)


gammas = ['0.04','0.12']
k24s = ['0.0','1.0']

for js,pair in enumerate(zip(gammas,k24s)):

    if js == 0:
        Lambdalist=['184']
        mtypes = ["s"]
        ltypes = ["--"]
    elif js == 1:
        Lambdalist=['27']
        mtypes = ["X","D"]
        ltypes = ["-.",":"]
        ms_lower = 43
        ms_upper = 12

    i_Lambda = int(Lambdalist[0])

    gamma = pair[0]
    k24 = pair[1]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24    
    scan['\\omega']=str(omega)

    obsfwd = ObservableData(["\\Lambda"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)

    Lambdas = obsfwd.data[:,0]

    x,y = obsfwd.surfacetwist()[1:],2*np.pi/obsfwd.eta()[1:]

    ax.plot(x,y,'.',color=colors[js])
    ax.plot(x[-1],y[-1],'rs')


dums = np.linspace(0.05,0.25,num=150,endpoint=True)
ax.plot(dums,np.cos(dums),'--')
    

ax.set_xlabel(r"$\psi(R)$",fontsize = 10)
ax.set_ylabel(r"$2\pi/\eta$",fontsize=10)


 


 
fig.subplots_adjust(left=0.1,right=0.9,bottom=0.15,top=0.95)
fig.savefig(obsfwd.observable_sname("twistvsdband",plot_format="pdf"))



plt.show()




import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from readparams import ReadParams

if len(sys.argv)<4:

    user_input = input("input string of a gamma,k24,omega values, "
                       "using comma as delimiter: ")
    gamma,k24,omega = user_input.split(',')

else:

    gamma,k24,omega = sys.argv[1],sys.argv[2],sys.argv[3]



width  = 3.487
height = width

#configure_fig_settings()

observable_list = ['E','R','eta','delta','surfacetwist']

fig = {}
ax = {}

data2d = {}

colors = sns.color_palette()

loadsuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]
savesuf = loadsuf



Lambdas = np.linspace(0,40,num=400,endpoint=True)


scan = {}
scan['\\gamma_s']=gamma
scan['k_{24}']=k24    
scan['\\omega']=str(omega)

obs = ObservableData(["\\Lambda"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                     savesuf=savesuf)


for observable in observable_list:


    fig[observable] = plt.figure()
    ax[observable] = fig[observable].add_subplot(1,1,1)

    fig[observable].set_size_inches(width,height)

    if observable == 'surfacetwist':
        ylabel = r'$\psi(R)$'
    elif observable == 'eta':
        ylabel = r'$2\pi/\eta$'
    elif observable == 'delta':
        ylabel = r'$\delta/\delta_0$'
    elif len(observable) > 1:
        ylabel = fr'$\{observable}$'
    else:
        ylabel = fr'${observable}$'


    ax[observable].plot(Lambdas,obs.R(),'.')

 
#for observable in observable_list:
    
#    fig[observable].subplots_adjust(left=0.2,right=0.95)
#    fig[observable].savefig(obs.observable_sname(observable))

plt.show()




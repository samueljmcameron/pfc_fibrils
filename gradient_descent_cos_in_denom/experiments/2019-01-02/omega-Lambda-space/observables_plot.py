import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../modules_gammak24/')
from plotobservables import PlotObservables
from gridmbyn import GridmByn
from readparams import ReadParams


width  = 3.487
height = width

#configure_fig_settings()

observable_list = ['E','R','eta','delta','surfacetwist']

observable_levels = {}

Lambda_yticks = {}

omega_yticks = {}

observable_levels['E'] = np.array([-4,-3,-2,-1,0],float)

observable_levels['R'] = np.array([0.225,0.23,1.9,2.0],float)

observable_levels['eta'] = np.array([0.985,0.99,0.993],float)

observable_levels['delta'] = np.array([0.7,0.75,0.8,0.81],float)

observable_levels['surfacetwist'] = np.array([0.16,0.18,0.198],float)

Lambda_yticks['surfacetwist'] = [0.2,0.3]

omega_yticks['surfacetwist'] = [0.2,0.3]

Lambda_yticks['delta'] = [0.805,0.81,0.815]

omega_yticks['delta'] = [0,0.4,0.8]

Lambda_yticks['eta'] = [0.985,0.99,0.995]

omega_yticks['eta'] = [0.993,0.994]

Lambda_yticks['R'] = [1,2]

omega_yticks['R'] = [0.25,0.35]

Lambda_yticks['E'] = [-2.5,-2]

omega_yticks['E'] = [-4,-2,0]


fig = {}
grid = {}

data2d = {}

colors = sns.color_palette()

savesuf = ["K_{33}","d_0"]
loadsuf = ["K_{33}","k_{24}","d_0","\\omega","\\gamma_s"]


num_Lambdas = 1000
num_omegas = 21
Lambdas = np.linspace(1,1000,num=num_Lambdas,endpoint=True)
omegas = np.linspace(10,30,num=num_omegas,endpoint=True)

for observable in observable_list:
    
    fig[observable] = plt.figure()

    fig[observable].set_size_inches(3*width,3*height)

    grid[observable] = GridmByn(fig[observable])


    data2d[observable]=np.empty([num_omegas,num_Lambdas],float)



for k24 in ['0.9','0.5','0.1']:

    for gamma in ['0.04','0.08','0.12']:

        print(f"gamma,k24 = {gamma},{k24}")
        scan = {}
        scan['\\gamma_s']=gamma
        scan['k_{24}']=k24



        # load in 2d grid of data in data2d for each observable at the
        # specified gamma,k24 pair.
        for om,omega in enumerate(omegas):

            scan['\\omega']=str(omega)

            rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

            obs = PlotObservables(["\\Lambda"],rp,scan_dir='scanforward')
            print(omega)
            for observable in observable_list:
        
                index = obs.ylabelstr_to_column(observable)
                data2d[observable][om,:] = obs.data[:,index]
                if observable == 'eta':
                    data2d[observable][om,:] = 2*np.pi/obs.data[:,index]
                    data2d[observable][om,0] = np.nan



        xlabel = r'$\Lambda$'
        ylabel = r'$\omega$'



        for observable in observable_list:

            if observable == 'surfacetwist':
                plabel = r'$\psi(R)$'
            elif observable == 'eta':
                plabel = r'$2\pi/\eta$'
            elif len(observable) > 1:
                plabel = fr'$\{observable}$'
            else:
                plabel = fr'${observable}$'



            z = data2d[observable]
            mask = np.zeros_like(z,dtype = bool)

            mask[:,0:2] = True
            mask[0:2,:] = True
            mask[:,24:25] = True
            mask[0:4,20:] = True
            mask[4:8,25:28] = True

            zmask = np.ma.array(z,mask=mask)

            cs = grid[observable].axarr[i][j]['main'].contourf(Lambdas,omegas,z)
            css = grid[observable].axarr[i][j]['main'].contour(Lambdas,omegas,zmask,
                                                               observable_levels[observable],
                                                               colors='r')
            grid[observable].axarr[i][j]['main'].clabel(css,fontsize=9,inline=1)

            cutout_omega_at_0 = data2d[observable][0,:]

            cutout_omega_at_15 = data2d[observable][15,:]

            cutout_omega_at_30 = data2d[observable][30,:]

            cutout_Lambda_at_0 = data2d[observable][:,0]

            cutout_Lambda_at_15 = data2d[observable][:,15]

            cutout_Lambda_at_30 = data2d[observable][:,30]


            grid[observable].axarr[i][j]['slice_const_y'].plot(Lambdas,cutout_omega_at_15,'.',
                                                               color='orange')
            #plt.setp(grid[observable].axarr[i][j]['slice_const_y'].get_xticklabels(),visible=False)
            grid[observable].axarr[i][j]['slice_const_y'].set_yticks(Lambda_yticks[observable])
            grid[observable].axarr[i][j]['slice_const_y'].set_ylabel(plabel)
            grid[observable].axarr[i][j]['slice_const_y'].set_xlim(Lambdas[0],Lambdas[-1])
            grid[observable].axarr[i][j]['slice_const_y'].set_xticks([0,5,10,15,20,25,30])

            ax_main[observable].get_xticklabels()[3].set_color("magenta")

            grid[observable].axarr[i][j]['slice_const_x'].plot(cutout_Lambda_at_15,omegas,'.',
                                      color='magenta')
            grid[observable].axarr[i][j]['slice_const_x'].set_xticks(omega_yticks[observable])
            grid[observable].axarr[i][j]['slice_const_x'].xaxis.tick_top()
            plt.setp(grid[observable].axarr[i][j]['slice_const_x'].get_yticklabels(),visible=False)
            grid[observable].axarr[i][j]['slice_const_x'].xaxis.set_label_position('top')
            grid[observable].axarr[i][j]['slice_const_x'].set_xlabel(plabel)
            grid[observable].axarr[i][j]['slice_const_x'].set_ylim(omegas[0],omegas[-1])

            grid[observable].axarr[i][j]['main'].get_yticklabels()[3].set_color("orange")

            grid[observable].axarr[i][j]['main'].set_xlabel(xlabel)
            grid[observable].axarr[i][j]['main'].set_ylabel(ylabel)

for observable in observable_list:
    
    fig[observable].subplots_adjust(left=0.18,right=0.95)
    fig[observable].savefig(obs.observable_sname(observable))

plt.show()




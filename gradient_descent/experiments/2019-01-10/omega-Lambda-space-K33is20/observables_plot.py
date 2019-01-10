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

omega_ylims = {}

Lambda_xlims = {}

vmins = {}

vmaxs = {}

vmins['E'] = -5

vmins['R'] = 0.01

vmins['surfacetwist'] = 0.1

vmins['eta'] = 0.98

vmins['delta'] = 0.8

vmaxs['E'] = 0

vmaxs['R'] = 1.0

vmaxs['surfacetwist'] = 0.3

vmaxs['eta'] = 1.0

vmaxs['delta'] = 1.0


omega_ylims['E'] = Lambda_xlims['E'] = [-6,1]

omega_ylims['R'] = Lambda_xlims['R'] = [0.02,5]

omega_ylims['surfacetwist'] = Lambda_xlims['surfacetwist'] = [0,0.3]

omega_ylims['eta'] = Lambda_xlims['eta'] = [0.97,1.0]

omega_ylims['delta'] = Lambda_xlims['delta'] = [0.85,1.0]


observable_levels['E'] = np.array([-4,-3,-2,-1,0],float)

observable_levels['R'] = np.array([0.05,0.1,0.5,1.0],float)

observable_levels['eta'] = np.array([0.97,0.98,0.99,1.0],float)

observable_levels['delta'] = np.array([0.87,0.9,0.93,0.95,0.99],float)

observable_levels['surfacetwist'] = np.array([0.1,0.15,0.2,0.25,0.3],float)

Lambda_yticks['surfacetwist'] = [0.05,0.15,0.25]

omega_yticks['surfacetwist'] = [0.05,0.15,0.25]

Lambda_yticks['delta'] = [0.9,0.95,1.0]

omega_yticks['delta'] = [0.9,0.95,1.0]

Lambda_yticks['eta'] = [0.98,0.99,1.0]

omega_yticks['eta'] = [0.98,0.99,1.0]

Lambda_yticks['R'] = [0.1,1.0]

omega_yticks['R'] = [0.1,1.0]

Lambda_yticks['E'] = [-5,-1,3]

omega_yticks['E'] = [-5,-1,3]


fig = {}
grid = {}

data2d = {}

colors = sns.color_palette()

savesuf = ["K_{33}","d_0"]
loadsuf = ["K_{33}","k_{24}","d_0","\\omega","\\gamma_s"]


num_Lambdas = max_Lambda = 1000
num_omegas = 30
Lambdas = np.linspace(1,max_Lambda,num=num_Lambdas,endpoint=True)
omegas = np.linspace(1,30,num=num_omegas,endpoint=True)

for observable in observable_list:
    
    fig[observable] = plt.figure()

    fig[observable].set_size_inches(3*width,3*height)

    grid[observable] = GridmByn(fig[observable])
    grid[observable].build_subplots_subgrid2x2(width_ratios=[3,1.2],
                                               hspace=0.1,wspace=0.1)

    data2d[observable]=np.empty([num_omegas,num_Lambdas],float)



for i,k24 in enumerate(['0.9','0.5','0.1']):

    for j,gamma in enumerate(['0.04','0.08','0.12']):

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
                datalength = len(obs.data[:,index])
                if datalength < num_Lambdas-1:
                    dumarray = np.ones(num_Lambdas-datalength)*np.nan
                    dataarray = np.concatenate((obs.data[:,index],dumarray))
                else:
                    dataarray = obs.data[:num_Lambdas,index]
                    
                if observable == 'eta':
                    data2d[observable][om,:] = 2*np.pi/dataarray
                    if Lambdas[0] == 0.0:
                        data2d[observable][om,0] = np.nan
                elif observable == 'delta':
                    data2d[observable][om,:] = dataarray/np.sqrt(2/3)
                else:
                    data2d[observable][om,:] = dataarray



        xlabel = r'$\Lambda$'
        ylabel = r'$\omega$'



        for observable in observable_list:

            if observable == 'surfacetwist':
                plabel = r'$\psi(R)$'
            elif observable == 'eta':
                plabel = r'$2\pi/\eta$'
            elif observable == 'delta':
                plabel = r'$\delta/\delta_0$'
            elif len(observable) > 1:
                plabel = fr'$\{observable}$'
            else:
                plabel = fr'${observable}$'



            z = data2d[observable]
            energies = data2d['E']
            mask = np.zeros_like(z,dtype = bool)

            """
            mask[:,0:2] = True
            mask[0:2,:] = True
            mask[:,24:25] = True
            mask[0:4,20:] = True
            mask[4:8,25:28] = True
            """
            zmask = z#np.ma.array(z,mask=mask)

            cs = grid[observable].axarr[i][j]['main'].contourf(Lambdas,omegas,z,
                                                               vmin = vmins[observable],
                                                               vmax = vmaxs[observable])
            css = grid[observable].axarr[i][j]['main'].contour(Lambdas,omegas,zmask,
                                                               observable_levels[observable],
                                                               colors='r')
            cssE = grid[observable].axarr[i][j]['main'].contour(Lambdas,omegas,energies,
                                                               [0],
                                                               colors='w')
            grid[observable].axarr[i][j]['main'].clabel(cssE,cssE.levels,inline=True,
                                                        manual = [(100,3)],
                                                        fmt={cssE.levels[0]:'E=0'})
            grid[observable].axarr[i][j]['main'].clabel(css,fontsize=9,inline=1)

            cutout_omega_at_15 = data2d[observable][15,:]

            cutout_Lambda_at_10 = data2d[observable][:,10]


            grid[observable].axarr[i][j]['slice_const_y'].plot(Lambdas,cutout_omega_at_15,'.',
                                                               color='orange')
            
            grid[observable].axarr[i][j]['slice_const_y'].set_yticks(Lambda_yticks[observable])

            grid[observable].axarr[i][j]['slice_const_y'].set_xlim(Lambdas[0],Lambdas[-1])
            grid[observable].axarr[i][j]['slice_const_y'].set_ylim(*omega_ylims[observable])
            if observable == 'R':
                grid[observable].axarr[i][j]['slice_const_y'].set_yscale('log')
                #labels = [tick.get_text() for tick in
                #          grid[observable].axarr[i][j]['slice_const_y'].get_yticklabels()]


            grid[observable].axarr[i][j]['main'].set_xscale('log')
            grid[observable].axarr[i][j]['main'].set_xticks([1,10,100,1000])
            grid[observable].axarr[i][j]['main'].get_xticklabels()[1].set_color("magenta")

            grid[observable].axarr[i][j]['slice_const_x'].plot(cutout_Lambda_at_10,omegas,'.',
                                      color='magenta')
            grid[observable].axarr[i][j]['slice_const_x'].set_xticks(omega_yticks[observable])
            grid[observable].axarr[i][j]['slice_const_x'].xaxis.tick_top()
            plt.setp(grid[observable].axarr[i][j]['slice_const_x'].get_yticklabels(),visible=False)
            grid[observable].axarr[i][j]['slice_const_x'].xaxis.set_label_position('top')

            grid[observable].axarr[i][j]['slice_const_x'].set_xlim(*Lambda_xlims[observable])
            grid[observable].axarr[i][j]['slice_const_x'].set_ylim(omegas[0],omegas[-1])
            if observable == 'R':
                grid[observable].axarr[i][j]['slice_const_x'].set_xscale('log')

            grid[observable].axarr[i][j]['slice_const_x'].set_xticklabels([])


            grid[observable].axarr[i][j]['main'].get_yticklabels()[3].set_color("orange")

            if i == 0:
                grid[observable].axarr[i][j]['slice_const_x'].set_xlabel(plabel)
            elif i == 2:
                grid[observable].axarr[i][j]['main'].set_xlabel(xlabel)
                grid[observable].axarr[i][j]['main'].annotate(rf'$\gamma={gamma}$',
                                                              xy=(0.7,-0.39),
                                                              xytext=(0.7,-0.4),
                                                              xycoords='axes fraction', 
                                                              ha='center',
                                                              va='center',
                                                              bbox=dict(boxstyle='square',
                                                                        fc='white'),
                                                              arrowprops=dict(arrowstyle='-[, widthB=9.5, lengthB=1.0', lw=2.0))
                                                              
            if j == 0:
                grid[observable].axarr[i][j]['main'].set_ylabel(ylabel)
                grid[observable].axarr[i][j]['slice_const_y'].set_ylabel(plabel)
                grid[observable].axarr[i][j]['main'].annotate(rf'$k_{{24}}={k24}$',
                                                              xy=(-0.49,0.7),
                                                              xytext=(-0.5,0.7),
                                                              xycoords='axes fraction', 
                                                              ha='center',
                                                              va='center',
                                                              bbox=dict(boxstyle='square',
                                                                        fc='white'),
                                                              arrowprops=dict(arrowstyle='-[, widthB=9.5, lengthB=1.0', lw=2.0))

            else:
                plt.setp(grid[observable].axarr[i][j]['slice_const_y'].get_yticklabels(),
                         visible=False)

for observable in observable_list:
    
    fig[observable].subplots_adjust(left=0.2,right=0.95)
    fig[observable].savefig(obs.observable_sname(observable))

plt.show()




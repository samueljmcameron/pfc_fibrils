#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../scripts')
from var_scan2d import loadfile_list, load_const_params
from universe_funcs import latex2string,string2latex
from variable_positions import return_position
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

def print_vals_near_min(arr,x_index,y_index):
    print("\t\t%10.5e"%arr[y_index+1,x_index])
    print("%10.5e\t%10.5e\t%10.5e"%(arr[y_index,x_index-1],arr[y_index,x_index],
                             arr[y_index,x_index+1]))
    print("\t\t%10.5e"%arr[y_index-1,x_index])

    return

def check_data(str_,string2compare2):
    if str_ != string2compare2:
        print("const_params.txt does not match the data"
              " that I think I'm plotting.")
        print("Expected %s, but got %s."%(string2compare2,str_))
        exit(1)
    return

def load_plt_array(str_,load_p):

    edited_str_ = latex2string(str_)
    var_array = np.loadtxt(load_p + '%ss.dat'%edited_str_)
    if var_array.size > 1:
        print(var_array)
        var_array = np.sort(var_array)
        np.savetxt(load_p + '%ss.dat'%edited_str_,var_array,
                   fmt = '%1.4e')
    else:
        var_array = np.array([var_array],float)
    return var_array


def ax_config(xlabel,ylabel,xscale,yscale,ax):

    # Set the basic axis properties. Mainly just made this
    # function for neatness.

    ax.set_ylabel(r"$%s$"%ylabel)
    ax.set_xlabel(r"$%s$"%xlabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.legend(frameon=False)

    return

def mk_mesh(d,params,excluded_param_position,numx,numy):

    xname = string2latex(d['scan_what_x'])
    yname = string2latex(d['scan_what_y'])

    x_position = return_position(xname)
    y_position = return_position(yname)

    if x_position > excluded_param_position:
        x_position -= 1

    if y_position > excluded_param_position:
        y_position -= 1

    x0 = params[x_position]
    xf = d['upperbound_x']
    print(numx)
    xs = np.linspace(x0,xf,num=numx,endpoint=True)

    y0 = params[y_position]
    yf = d['upperbound_y']+0.001
    ys = np.linspace(y0,yf,num=numy,endpoint=True)
    print(numy)
    return np.meshgrid(xs,ys)
    

def plot_scanE2d(fig,ax,colors,d,var,params,var_position,
                 varied_param_name,load_p):

    # plots E vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.


    load_str = loadfile_list(params,var,var_position)
    fname= ("%s_%s_%s_energy_%s")%(load_p,d['scan_what_x'],
                                   d['scan_what_y'],load_str)
    
    legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                      var)
    
    # load and plot E vs x
    EE = np.loadtxt(fname + ".txt")

    numy,numx = EE.shape

    xx,yy = mk_mesh(d,params,var_position,numx,numy)

    cs = ax.imshow(EE,vmin=EE.min(),vmax=EE.max(),
                   cmap='Reds_r',#norm=LogNorm(),
                   extent=[xx.min(),xx.max(),yy.min(),yy.max()],
                   aspect=xx.max()/yy.max(),origin='lower')

    cbar = fig.colorbar(cs,ax=ax,fraction=0.046,pad=0.04)#cax=cax)
    cs.set_clim(EE.min(),EE.max())

    y_index,x_index = np.unravel_index(EE.argmin(),EE.shape)

    ax.plot(xx[x_index,y_index],yy[x_index,y_index],'ko')

    print(x_index,y_index)

    print_vals_near_min(EE,x_index,y_index)

    ax.set_xlim(xx.min(),xx.max())
    ax.set_ylim(yy.min(),yy.max())

    return cs

def plot_scanderivEx(fig,ax,colors,d,var,params,var_position,
                     varied_param_name,load_p,which_deriv):

    # plots dEdx vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.

    load_str = loadfile_list(params,var,var_position)
    fname= ("%s_%s_%s_deriv_"
            "energy_%s_%s")%(load_p,d['scan_what_x'],
                             d['scan_what_y'],which_deriv,
                             load_str)
    
    legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                      var)
    
    # load and plot E vs x
    dEEdxx = np.loadtxt(fname + ".txt")

    fname= ("%s_%s_%s_energy_%s")%(load_p,d['scan_what_x'],
                                   d['scan_what_y'],load_str)

    EE = np.loadtxt(fname + ".txt")

    numy,numx = dEEdxx.shape

    xx,yy = mk_mesh(d,params,var_position,numx,numy)

    cs = ax.imshow(dEEdxx,vmin=-1,vmax=1,#dEEdxx.min(),vmax=dEEdxx.max(),
                   cmap=cm.cool,
                   extent=[xx.min(),xx.max(),yy.min(),yy.max()],
                   aspect=xx.max()/yy.max(),origin='lower')

    cbar = fig.colorbar(cs,ax=ax,fraction=0.046,pad=0.04)
    cs.set_clim(-1,1)

    y_index,x_index = np.unravel_index(EE.argmin(),EE.shape)
    
    print(x_index,y_index)

    ax.plot(xx[x_index,y_index],yy[x_index,y_index],'ko')

    print_vals_near_min(dEEdxx,x_index,y_index)

                

    return cs

def plot_scanpsi(ax,colors,d,var,params,
                 var_position,varied_param_name,load_p):

    # plots all psi vs r which minimize the energy,
    # for all the values of var in the array var_array. 
    # varied_param_name is the name (a string) of var. 
    # load_p is the path which the data for the plots is
    # loaded from. If there are multiple psi(r) which
    # minimize the energy, it is stated in the STDOUT.


    # length is the number of data points in the psi(r)
    # data, assuming only one unique configuration
    # minimizes the energy.

    length = 2*2*2*2*2*2*2*2*2*2*2+1


    load_str = loadfile_list(params,var,var_position)
    fname= ("%s_%s_%s_psivsr_%s")%(load_p,d['scan_what_x'],
                                   d['scan_what_y'],load_str)

    if (varied_param_name == 'gamma_s'
        or varied_param_name == 'Lambda'
        or varied_param_name == 'eta'):
        legend_label = "\%s=\SI{%1.1e}{}"%(varied_param_name,
                                           var)
    else:
        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                           var)


    # load and plot psi vs r
    psivsr = np.loadtxt(fname + ".txt")
    if psivsr.size != 0:
        print(psivsr.size)
        rs = psivsr[:,0]
        psis = psivsr[:,1]
        dpsidrs = psivsr[:,2]


        if len(rs) % length != 0:
            print("The length of psi(r) is not what "
                  "was expected!")
            print("Exiting to system while trying to plot!")
            exit(1)

        j = 1
        while (len(rs)-j*length>=0):
            if j == 1:
                ax.plot(rs[(j-1)*length:j*length],
                        psis[(j-1)*length:j*length],
                        '-',lw = 2,color = colors[i],
                        label = r"$%s$"%legend_label)
            else:
                ax.plot(rs[(j-1)*length:j*length],
                        psis[(j-1)*length:j*length],
                        '-',lw = 2,color = colors[i])
            if j > 1:
                print("multiple psi(r) configurations "
                      " with the same energy!")
            j += 1

    return

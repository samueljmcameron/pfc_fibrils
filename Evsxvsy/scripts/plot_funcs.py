#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../scripts')
from var_scan2d import loadfile_list, load_const_params
from universe_funcs import latex2string,string2latex
from variable_positions import return_position
from matplotlib import cm, ticker

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
        print(var_array.size)
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

def mk_mesh(d,params,excluded_param_position):

    var_x_spacing = 0.001
    var_y_spacing = 0.001

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
    numx = int(round((xf-x0)/var_x_spacing))
    print(numx)
    xs = np.linspace(x0,xf,num=numx,endpoint=True)

    y0 = params[y_position]
    yf = d['upperbound_y']
    numy = int(round((yf-y0)/var_y_spacing))
    ys = np.linspace(y0,yf,num=numy,endpoint=True)
    print(numy)
    return np.meshgrid(xs,ys)
    

def plot_scanE2d(ax,colors,d,var,params,var_position,
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

    print(EE.shape)

    xx,yy = mk_mesh(d,params,var_position)

    print(xx.shape)
    print(yy.shape)

    cs = ax.imshow(EE.T,vmin=EE.min(),vmax=EE.max(),
                   cmap=cm.coolwarm,
                   extent=[xx.min(),xx.max(),yy.min(),yy.max()],
                   aspect=1.0/5.0,origin='lower')

    return cs

def plot_scanderivEx(ax,colors,d,var,params,var_position,
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

    print(dEEdxx.shape)

    xx,yy = mk_mesh(d,params,var_position)

    print(xx.shape)
    print(yy.shape)

    cs = ax.imshow(dEEdxx.T,vmin=-10,vmax=10,#dEEdxx.min(),vmax=dEEdxx.max(),
                   cmap=cm.coolwarm,
                   extent=[xx.min(),xx.max(),yy.min(),yy.max()],
                   aspect=1.0/5.0,origin='lower')


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

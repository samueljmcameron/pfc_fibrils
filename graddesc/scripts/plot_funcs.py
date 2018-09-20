#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../scripts')
from var_scan import loadfile_list, load_const_params,latex2string

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


def plot_desc(what_plot,ax,colors,d,var_array,params,
              var_position,varied_param_name,load_p):

    # plots E vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params[:-3],var,var_position)
        fname= ("%s_%s_%s")%(load_p,what_plot,load_str)

        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                          var)

        # load and plot y vs x
        yvsx = np.loadtxt(fname + ".txt")
        xs = yvsx[:,0]
        ys = yvsx[:,1]

        ax.plot(xs,ys,'-o',lw = 2,color = colors[i],
                label = r"$%s$"%legend_label)

    return

def plot_descpsi(ax,colors,d,var_array,params,
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

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params[:-3],var,var_position)
        fname= ("%s_psivsr_%s")%(load_p,load_str)

        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                          var)

        # load and plot psi vs r
        psivsr = np.loadtxt(fname + ".txt")
        if psivsr.size != 0:
            print(psivsr.size)
            rs = psivsr[:,0]
            psis = psivsr[:,1]
            dpsidrs = psivsr[:,2]

            zero_index = np.argmin(rs[1:])+1
            if rs[zero_index] == 0:
                length = len(rs[:zero_index])
            else:
                length = len(rs)

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

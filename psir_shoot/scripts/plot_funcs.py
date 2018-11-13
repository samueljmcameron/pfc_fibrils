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

def load_plt_array(str_,load_p,suffix=''):

    edited_str_ = latex2string(str_)
    var_array = np.loadtxt(load_p + '%ss%s.dat'%(edited_str_,suffix))
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


def plot_psivsr(ax,colors,d,var_array,params,
                var_position,varied_param_name,load_p):

    # plots E vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.

    max_psip0s = np.empty([0],float)

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params,var,var_position)
        fname= ("%s_%s_%s")%(load_p,"psivsr",load_str)

        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                          var)

        # load and plot y vs x
        data = np.loadtxt(fname + ".txt")
        rs = data[:,0]
        psis = data[:,1]
        
        ax.plot(rs,psis,'-',lw = 2,color = colors[i],
                label = r"$%s$"%legend_label)

    return

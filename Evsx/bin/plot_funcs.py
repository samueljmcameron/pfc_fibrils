#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../bin')
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


def plot_scanE(ax,colors,d,var_array,params,
               var_position,varied_param_name,load_p):

    # plots E vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params,var,var_position)
        fname= ("%s_%s_energy_%s")%(load_p,d['scan_what'],
                                    load_str)


        #if (varied_param_name == 'Lambda'
        #    or varied_param_name == 'omega'
        #    or varied_param_name == 'eta'
        #    or varied_param_name == 'delta'
        #    or varied_param_name == 'gamma_s'
        #    or varied_param_name == 'gamma_t'):
        #    legend_label = "\%s=\SI{%1.1e}{}"%(varied_param_name,
        #                                       var)
        #        else:
        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                               var)

        # load and plot E vs x
        Evsx = np.loadtxt(fname + ".txt")
        xs = Evsx[:,0]
        Es = Evsx[:,1]
        dEdxs = Evsx[:,2]
        psixs = Evsx[:,3]

        ax.plot(xs,Es,'-',lw = 2,color = colors[i],
                label = r"$%s$"%legend_label)
        for i,x in enumerate(xs[1:]):
            if dEdxs[i]*dEdxs[i-1] <= 0 and dEdxs[i] >0:
                ax.plot(xs[i],Es[i],'k.')

    return

def plot_scanderivE(ax,colors,d,var_array,params,
                    var_position,varied_param_name,load_p):

    # plots E vs d['scan_what'], for all the values of 
    # var in the array var_array. varied_param_name is
    # the name (a string) of var. load_p is the path
    # which the data for the plots is loaded from.

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params,var,var_position)
        fname= ("%s_%s_energy_%s")%(load_p,d['scan_what'],
                                    load_str)


        #if (varied_param_name == 'Lambda'
        #    or varied_param_name == 'omega'
        #    or varied_param_name == 'eta'
        #    or varied_param_name == 'delta'
        #    or varied_param_name == 'gamma_s'
        #    or varied_param_name == 'gamma_t'):
        #    legend_label = "\%s=\SI{%1.1e}{}"%(varied_param_name,
        #                                       var)
        #        else:
        legend_label = "%s=\SI{%1.1e}{}"%(varied_param_name,
                                               var)

        # load and plot E vs x
        Evsx = np.loadtxt(fname + ".txt")
        xs = Evsx[:,0]
        Es = Evsx[:,1]
        dEdxs = Evsx[:,2]
        psixs = Evsx[:,3]

        ax.plot(xs,dEdxs,'-',lw = 2,color = colors[i],
                label = r"$%s$"%legend_label)

    return


def plot_scanpsi(ax,colors,d,var_array,params,
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

    for i,var in enumerate(var_array):

        load_str = loadfile_list(params,var,var_position)
        fname= ("%s_%s_psivsr_%s")%(load_p,d['scan_what'],
                                    load_str)

        if (varied_param_name == 'gamma_s'
            or varied_param_name == 'gamma_t'
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

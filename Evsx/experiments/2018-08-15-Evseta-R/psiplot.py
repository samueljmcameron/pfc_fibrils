#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../../bin')
sys.path.append('../../bin')
from fig_settings import configure_fig_settings
from var_scan import loadfile_list, load_const_params,latex2string
from plot_funcs import plot_scanpsi, ax_config, check_data
from plot_funcs import load_plt_array
import seaborn as sns


if __name__=='__main__':

    configure_fig_settings() # keep

    d,params,str_,var_position = load_const_params()

    check_data(str_,'\gamma_t')

    save_p = "results/"
    load_p = "data/"

    var_array = load_plt_array(str_,load_p)

    fig,ax = plt.subplots() # keep
    width  = 3.487 # keep
    height = width # keep
    fig.set_size_inches(width,height) # keep

    colors = sns.color_palette('muted',len(var_array)) # keep

    plot_scanpsi(ax,colors,d,var_array,params,var_position,
                 str_,load_p)

    ylabel = "\psi(r)" # keep
    xlabel = 'r' # keep
    xscale = 'linear' # keep
    yscale = 'linear' # keep

    ax_config(xlabel,ylabel,xscale,yscale,ax) # keep

    edited_str_ = latex2string(str_)

    sname = "psi_%svsr-%ss"%(d['scan_what'],
                             edited_str_) # keep

    fig.tight_layout() # keep
    fig.savefig(save_p + sname + ".pdf") # keep

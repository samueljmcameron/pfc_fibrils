#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../../scripts')
sys.path.append('../../scripts')
from fig_settings import configure_fig_settings
from var_scan2d import loadfile_list, load_const_params
from universe_funcs import latex2string, string2latex
from plot_funcs import plot_scanpsi, ax_config, check_data
from plot_funcs import load_plt_array
import seaborn as sns


if __name__=='__main__':

    configure_fig_settings() # keep

    d,params,str_,var_position = load_const_params()

    check_data(str_,sys.argv[1])

    index_var = int(sys.argv[2])

    save_p = "results/"
    load_p = "data/"

    var_array = load_plt_array(str_,load_p)

    fig,ax = plt.subplots() # keep
    width  = 3.487 # keep
    height = width # keep
    fig.set_size_inches(width,height) # keep

    colors = sns.color_palette('muted',len(var_array)) # keep

    plot_scanpsi(ax,colors,d,var_array[index_var],params,var_position,
                 str_,load_p)

    ylabel = "\psi(r)" # keep
    xlabel = 'r' # keep
    xscale = 'linear' # keep
    yscale = 'linear' # keep

    ax_config(xlabel,ylabel,xscale,yscale,ax) # keep

    edited_str_ = latex2string(str_)

    sname = "psi_%svs%svsr-%ss-%1.2e"%(d['scan_what_x'],d['scan_what_y'],
                                       edited_str_,var_array[index_var])

    fig.tight_layout() # keep
    fig.savefig(save_p + sname + ".pdf") # keep
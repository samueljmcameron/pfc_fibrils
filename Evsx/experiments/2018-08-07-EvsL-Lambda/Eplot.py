#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../../bin')
sys.path.append('../../bin')
from fig_settings import configure_fig_settings
from var_scan import loadfile_list, load_const_params
from plot_funcs import plot_scanE, ax_config, check_data
from plot_funcs import load_plt_array
import seaborn.apionly as sns


if __name__=='__main__':

    configure_fig_settings()


    d,params,str_,var_position = load_const_params()

    check_data(str_,'Lambda')


    save_p = "results/"
    load_p = "data/"

    var_array = load_plt_array(str_,load_p)

    colors = sns.color_palette('muted',len(var_array))

    fig,ax = plt.subplots()
    width  = 3.487
    height = width
    fig.set_size_inches(width,height)

    plot_scanE(ax,colors,d,var_array,params,var_position,
               str_,load_p)

    ylabel = "E"
    xlabel = d['scan_what']
    xscale = 'linear'
    yscale = 'log'

    ax_config(xlabel,ylabel,xscale,yscale,ax)


    sname = "Evs%s-%ss"%(d['scan_what'],
                        str_)

    fig.tight_layout()
    fig.savefig(save_p + sname + ".pdf")

#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../../scripts')
sys.path.append('../../scripts')
from fig_settings import configure_fig_settings
from var_scan import loadfile_list, load_const_params,latex2string
from plot_funcs import plot_obs, ax_config, check_data
from plot_funcs import load_plt_array
import seaborn as sns


if __name__=='__main__':

    configure_fig_settings()


    d,params,str_,var_position = load_const_params()

    check_data(str_,sys.argv[1])


    save_p = "results/"
    load_p = "data/"

    var_array = load_plt_array(str_,load_p)

    obss = ['energy','\delta','\eta','R','surfacetwist']

    colors = sns.color_palette('muted',len(obss))


    for i,ylabel in enumerate(obss):
        edited_ylabel = latex2string(ylabel)

        fig,ax = plt.subplots()
        width  = 3.487
        height = width
        fig.set_size_inches(width,height)

        if ylabel != 'energy' and ylabel != 'surfacetwist':
            plot_obs('dEd%s'%edited_ylabel,ax,colors[i],d,
                     var_array,params,var_position,str_,load_p)
        else:
            plot_obs(edited_ylabel,ax,colors[i],d,var_array,
                     params,var_position,str_,load_p)

        xlabel = sys.argv[1]
        xscale = 'log'
        yscale = 'linear'

        if ylabel == 'surfacetwist':
            ylabel = '\psi(R)'
        elif ylabel == 'energy':
            ylabel = 'E'

        ax_config(xlabel,ylabel,xscale,yscale,ax)

        edited_str_ = latex2string(str_)

        sname = "%svs%ss"%(edited_ylabel,edited_str_)

        fig.tight_layout()
        fig.savefig(save_p + sname + ".pdf")

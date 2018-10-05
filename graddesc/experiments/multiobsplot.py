#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
sys.path.append('../../../scripts')
sys.path.append('../../scripts')
from variable_positions import return_position
from fig_settings import configure_fig_settings
from var_scan import loadfile_list, load_const_params,latex2string
from plot_funcs import plot_obs, ax_config, check_data
from plot_funcs import load_plt_array
import seaborn as sns


if __name__=='__main__':

    configure_fig_settings()

    nmbr_of_second_variable_values = int(sys.argv[3])
    listofdates = sys.argv[4].split(',')
    print(listofdates)

    colors = sns.color_palette('muted',nmbr_of_second_variable_values)

    for ylabel in ['energy','\delta','\eta','R','surfacetwist']:

        fig,ax = plt.subplots()
        width  = 3.487
        height = width
        fig.set_size_inches(width,height)

        edited_ylabel = latex2string(ylabel)

        for i in range(nmbr_of_second_variable_values):
            d,params,str_,var_position = load_const_params(i)
            edited_str_ = latex2string(str_)
            check_data(str_,sys.argv[1])
                    
            save_p = "results/"
            load_p = "../%s-%s-%sis%.3lf/data/"%(listofdates[i],edited_str_,
                                                 latex2string(sys.argv[2]),
                                                 d[sys.argv[2]])
        
            var_array = load_plt_array(str_,load_p)
            
            
            label = r'$%s=\SI{%1.1e}{}$'%(sys.argv[2],d[sys.argv[2]])

            if edited_ylabel != 'energy' and edited_ylabel != 'surfacetwist':
                plot_obs('dEd%s'%edited_ylabel,ax,colors[i],d,var_array,params,var_position,
                         str_,load_p,label=label)
            else:
                if ylabel == 'surfacetwist':
                    ylabel = '\psi(R)'
                elif ylabel == 'energy':
                    ylabel = 'E'
                plot_obs(edited_ylabel,ax,colors[i],d,var_array,params,var_position,
                         str_,load_p,label=label)            

        xlabel = sys.argv[1]
        xscale = 'log'
        yscale = 'linear'

        ax_config(xlabel,ylabel,xscale,yscale,ax)
        ax.legend(frameon=False)

        edited_str_ = latex2string(str_)

        sname = "%svs%ss"%(edited_ylabel,edited_str_)

        fig.tight_layout()
        fig.savefig(save_p + sname + ".pdf")

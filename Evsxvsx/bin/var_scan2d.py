#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
import re
sys.path.append('../../../bin')
from variable_positions import return_position
from universe_funcs import not_number, load_const_params
from universe_funcs import loadfile_list,latex2string


def argv_list(init_path,params,var,var_position,
              scan_what_x,scan_what_y):

    # Load in the list of parameters which must be fed to
    # the executable scan (with src code from scan.c).

    args = ("%s %s %lf %s "
            "%s %s")%(init_path,
                      ' '.join(map(str,params[0:var_position])),
                      var,
                      ' '.join(map(str,params[var_position:])),
                      scan_what_x,scan_what_y)
    return args

def one_run(var):

    # Compute E vs x for one value of varied_param.

    
    d,params,varied_param_name,var_position = load_const_params()

    successful_var = []
    init_path = "../../tmp_data/"
    
        
    # run scan executable for var
    variable_scanE2d(init_path,params,var,var_position,
                     d['scan_what_x'],d['scan_what_y'],
                     successful_var)

    # store whichever calculations were successful into
    # a .dat file in the results folder
    successful_var = np.array(successful_var,
                                   'float')

    if successful_var.size != 0:
        edited_name = latex2string(varied_param_name)
        with open('data/%ss.dat'%edited_name,'a') as file:
            file.write('{:1.4e}'.format(successful_var[0])
                       +'\n')

    return


def variable_scanE2d(init_path,params,var,var_position,
                     scan_what_x,scan_what_y,
                     successful_calc_list):

    # For a specified set of constant params, as well
    # as the parameter you are varying to see its
    # effect on 'Evs%s'%scan_what (i.e. var in the
    # argument list above), attempt to compute 
    # the data for E vs '%s'%scan_what. If it fails,
    # return the stdout of the run.


    args = argv_list(init_path,params,var,var_position,
                     scan_what_x,scan_what_y)

    cmd = "../../bin/scan " + args

    if(subprocess.call(cmd,shell=True,
                       stderr=subprocess.STDOUT)==0):
        successful_calc_list.append(var)
        load_str = loadfile_list(params,var,var_position)
        E_f = mk_fname(init_path,scan_what_x,scan_what_y,
                       'energy',load_str)
        psi_f = mk_fname(init_path,scan_what_x,scan_what_y,
                         'psivsr',load_str)
        dEdx_f = mk_fname(init_path,scan_what_x,scan_what_y,
                          'deriv_energy_%s'%scan_what_x,load_str)
        dEdy_f = mk_fname(init_path,scan_what_x,scan_what_y,
                          'deriv_energy_%s'%scan_what_y,load_str)
        srftwst_f = mk_fname(init_path,scan_what_x,scan_what_y,
                             'surfacetwist',load_str)


        cmdE = "mv " + E_f + " data/"
        subprocess.call(cmdE,shell=True)
        cmdpsi = "mv " + psi_f + " data/"
        subprocess.call(cmdpsi,shell=True)
        cmddEdx = "mv " + dEdx_f + " data/"
        subprocess.call(cmddEdx,shell=True)
        cmddEdy = "mv " + dEdy_f + " data/"
        subprocess.call(cmddEdy,shell=True)
        cmdsrftwst = "mv " + srftwst_f + " data/"
        subprocess.call(cmdsrftwst,shell=True)
    else:
        print(subprocess.STDOUT)
    return 

def mk_fname(init_path,scan_what_x,scan_what_y,
             m_str,load_str):
    
    fname = ("%s_%s_%s_%s_%s"
             ".txt")%(init_path,scan_what_x,
                      scan_what_y,m_str,load_str)
    return fname

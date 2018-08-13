#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
sys.path.append('../../../bin')
from variable_positions import return_position

def not_number(s):

    # Check if s is not in a format which can be translated
    # into a float, and if it is NOT, return True. Otherwise,
    # return False.

    try:
        float(s)
        return False
    except ValueError:
        return True

def load_const_params():

    # Load in dictionary of the constant parameter values,
    # put those values in the pre-specified order (in the 
    # 'const_params.txt' file), and then write the name of
    # the parameter which is being varied (e.g. 'gamma_t')
    # to the string variable varied_key_name. Also returns
    # the position of the variable which is being varied,
    # according to the function return_position(key).

    d = {}
    params = []
    varied_param_name = ""

    with open("data/constant_params.txt") as f:
        for line in f:
            (key,val) = line.split()
            if val == "#varying":
                d[key] = val
                varied_param_name = key
                var_position = return_position(key)
            elif not_number(val):
                d[key] = val
            else:
                d[key] = float(val)
                params.append(float(val))

    params = np.array(params,float)

    return d, params, varied_param_name, var_position

def argv_list(init_path,params,var,var_position,
              scan_what):

    # Load in the list of parameters which must be fed to
    # the executable scan (with src code from scan.c).

    args = ("%s %s %lf "
            "%s %s")%(init_path,
                      ' '.join(map(str,params[0:var_position])),
                      var,
                      ' '.join(map(str,params[var_position:])),
                      scan_what)
    return args

def loadfile_list(params,var,var_position):
    
    # Replicate the specific format of the file names which
    # are saved by the executable scan, but only the last
    # part (e.g. the part after 'energy_' or 'psivsr').

    s = "%s_%1.4e_%s"%('_'.join(map('{:1.4e}'.format,
                                    params[:var_position])),
                       var,
                       '_'.join(map('{:1.4e}'.format,
                                    params[var_position:])))
    return s

def one_run(var):

    # Compute E vs x for one value of varied_param.

    
    d,params,varied_param_name,var_position = load_const_params()

    successful_var = []
    init_path = "../../tmp_data/"
    
        
    # run scan executable for var
    variable_scanE(init_path,params,var,var_position,
                   d['scan_what'],successful_var)

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

def latex2string(str_):
    if str_[0] == '\\':
        edited = str_[1:]
    else:
        edited = str_
    return edited


def variable_scanE(init_path,params,var,var_position,
                   scan_what,successful_calc_list):

    # For a specified set of constant params, as well
    # as the parameter you are varying to see its
    # effect on 'Evs%s'%scan_what (i.e. var in the
    # argument list above), attempt to compute 
    # the data for E vs '%s'%scan_what. If it fails,
    # return the stdout of the run.


    args = argv_list(init_path,params,var,var_position,
                     scan_what)

    cmd = "../../bin/scan " + args

    if(subprocess.call(cmd,shell=True,
                       stderr=subprocess.STDOUT)==0):
        successful_calc_list.append(var)
        load_str = loadfile_list(params,var,var_position)
        Efile =  ("%s_%s_energy_%s"
                  ".txt")%(init_path,scan_what,
                           load_str)
        psifile = ("%s_%s_psivsr_%s"
                   ".txt")%(init_path,scan_what,
                            load_str)


        cmdE = "mv " + Efile + " data/"
        subprocess.call(cmdE,shell=True)
        cmdpsi = "mv " + psifile + " data/"
        subprocess.call(cmdpsi,shell=True)
    else:
        print(subprocess.STDOUT)
    return 


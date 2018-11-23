#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
import re
sys.path.append('../../../scripts')
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

def load_const_params(index=None):

    # Load in dictionary of the constant parameter values,
    # put those values in the pre-specified order (in the 
    # 'const_params.txt' file), and then write the name of
    # the parameter which is being varied (e.g. 'gamma_s')
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
            elif ',' in val:
                if index == None:
                    print("array of values for parameter %s, but "
                          "no index value to select which value "
                          "of the array is to be used."%key)
                    exit(1)
                array_val = np.array(val.split(','),float)[index]
                d[key] = array_val
                params.append(array_val)
            elif not_number(val):
                d[key] = val
            else:
                d[key] = float(val)
                params.append(float(val))

    params = np.array(params,float)

    return d, params, varied_param_name, var_position

def argv_list(init_path,params,var,var_position):

    # Load in the list of parameters which must be fed to
    # the executable graddesc (with src code from graddesc.c).

    args = ("%s %s %lf "
            "%s")%(init_path,
                      ' '.join(map(str,params[0:var_position])),
                      var,
                      ' '.join(map(str,params[var_position:])))
    return args

def loadfile_list(params,var,var_position):
    
    # Replicate the specific format of the file names which
    # are saved by the executable graddesc, but only the last
    # part (e.g. the part after 'energy_' or 'psivsr').

    if var_position == 0:
        s = "%1.4e_%s"%(var,
                        '_'.join(map('{:1.4e}'.format,
                                     params[var_position:])))
    else:
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
    
        
    # run shooting executable for var
    shooting(init_path,params,var,var_position,successful_var)

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
    edited = re.sub(r'\{',r'',re.sub(r'\}',r'',re.sub(r'\\',r'',str_)))
    return edited


def shooting(init_path,params,var,var_position,successful_calc_list):


    args = argv_list(init_path,params,var,var_position)

    cmd = "../../bin/psivsr " + args

    if(subprocess.call(cmd,shell=True,
                       stderr=subprocess.STDOUT)==0):

        successful_calc_list.append(var)
        load_str = loadfile_list(params[:-2],var,var_position)

        for fpart in ['psivsr']:
            file = ("%s_%s_%s"
                    ".txt")%(init_path,fpart,load_str)

            cmd = "mv " + file + " data/"
            #cmd = "mv " + file + " " + init_path + fpart + ".txt"
            subprocess.call(cmd,shell=True)

    else:
        print(subprocess.STDOUT)

    return 

#!/usr/bin/env python
from __future__ import print_function, division
import subprocess
import numpy as np
import sys
import re
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


def loadfile_list(params,var,var_position):
    
    # Replicate the specific format of the file names which
    # are saved by the executable scan, but only the last
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


def latex2string(str_):
    edited = re.sub(r'\{',r'',re.sub(r'\}',r'',re.sub(r'\\',r'',str_)))
    return edited

def string2latex(str_):
    if str_ == 'eta':
        return '\eta'
    elif str_ == 'delta':
        return '\delta'
    else:
        return str_

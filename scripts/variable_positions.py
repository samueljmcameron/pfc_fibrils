from __future__ import print_function, division
import numpy as np

def return_position(key):
    
    # Given a dictionary key, return the "position"
    # (as arbitrarily decided by me).

    if key == 'K_{33}':
        return 0
    elif key == 'k_{24}':
        return 1
    elif key == '\Lambda':
        return 2
    elif key == 'd_0':
        return 3
    elif key == '\omega':
        return 4
    elif key == 'R':
        return 5
    elif key == '\eta':
        return 6
    elif key == '\delta':
        return 7
    elif key == '\gamma_s':
        return 8
    else:
        print("no position could be found for"
              "the specified key!")
        exit(255)

    return

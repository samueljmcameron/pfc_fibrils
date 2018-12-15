from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import sys
from matplotlib import rc

def configure_fig_settings():

    font = {'family':'Times New Roman','weight':'normal','size':'9'}
    rc('font', **font)
    rc('text',usetex=True)
    params = {'text.latex.preamble':[r'\usepackage{siunitx}',
                                     r'\usepackage{amsmath}']}

    
    plt.rcParams.update(params)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

    return

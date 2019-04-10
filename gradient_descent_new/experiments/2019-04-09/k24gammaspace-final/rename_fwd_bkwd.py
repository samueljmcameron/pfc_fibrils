import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import os
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from gridmbyn import GridmByn
from readparams import ReadParams



configure_fig_settings()

observable_list = ['E','R','eta','delta','surfacetwist']

colors = sns.color_palette()

savesuf = ["K_{33}","\\Lambda","\\omega"]
loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]


num_k24s=51
k24s = np.linspace(-1,1,num=101,endpoint=True)


Lambda = sys.argv[1]
omega = sys.argv[2]

scan = {}
scan['\\Lambda']=Lambda
scan['\\omega']=omega

# load in 2d grid of data in data2d for each observable at the
# specified gamma,k24 pair.
for i,k24 in enumerate(k24s):

    scan['k_{24}']=str(k24)

    obsfwd = ObservableData(["\\gamma_s"],scan_dir='scanforward',scan=scan,
                            loadsuf=loadsuf,savesuf=savesuf)
    obsbkwd = ObservableData(["\\gamma_s"],scan_dir='scanbackward',scan=scan,
                             loadsuf=loadsuf,savesuf=savesuf)
    
    bkwdexists = os.path.isfile(obsbkwd.observables_fname())

    if not bkwdexists:

        # save a copy of scanforward data under the scanbackward filename,
        # as if the file doesn't exist it's assumed that the scanbackward
        # and scanforward files are the same (below the critical point region
        # in gamma,k24 space)

        np.savetxt(obsbkwd.observables_fname(),obsfwd.data,
                   fmt='\t'.join(["%13.6e"] + ["%15.8e"]*5))
    


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from gridmbyn import GridmByn
from readparams import ReadParams


width  = 3.487
height = width

#configure_fig_settings()

observable_list = ['E','R','eta','delta','surfacetwist']

colors = sns.color_palette()

savesuf = ["K_{33}","\\Lambda","\\omega"]
loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]


num_gammas=101
gammas = np.linspace(0.01,0.4,num=101,endpoint=True)
num_k24s=51
k24s = np.linspace(0,1,num=51,endpoint=True)





scan = {}
scan['\\Lambda']=sys.argv[1]
scan['\\omega']=sys.argv[2]



# load in 2d grid of data in data2d for each observable at the
# specified gamma,k24 pair.
for k24 in k24s:
 

    scan['k_{24}']=str(k24)

    obsfwd = ObservableData(["\\gamma_s"],scan_dir='scanforward',scan=scan,
                            loadsuf=loadsuf,savesuf=savesuf)
    obsbkwd = ObservableData(["\\gamma_s"],scan_dir='scanbackward',scan=scan,
                             loadsuf=loadsuf,savesuf=savesuf)

    obsbkwd.sort_observables()

    end_bkwd = obsbkwd.data.shape[0]

    datanew = np.concatenate((obsbkwd.data,obsfwd.data[end_bkwd:,:]),axis=0)

    np.savetxt(obsbkwd.observables_fname(),datanew,
               fmt='\t'.join(["%13.6e"]+["%15.8e"]*5))

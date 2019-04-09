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



configure_fig_settings()

observable_list = ['E','R','eta','delta','surfacetwist']

colors = sns.color_palette()

savesuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]
loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega"]


num_gammas=101
gammas = np.linspace(0.01,0.4,num=101,endpoint=True)
num_k24s=51

k24_0s = np.linspace(-1,0,num=50,endpoint=False)

k24s = np.linspace(0,1,num=51,endpoint=True)


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


    print(k24)
    # find all spots where two distinct phases exist (two different radii values)
    bool_1 = np.isclose(obsfwd.R(),obsbkwd.R(),rtol=1e-3)
    find_j = np.where(bool_1==False)[0]

    if find_j.size > 0:   # if two distinct phases exist, then:

        # smallest gamma value that they both exist at
        jsmall = find_j[0]

        # largest gamma value that they both exist at
        jlarge = find_j[-1]

        # find the point where the fwd E becomes larger than the bkwd E
        j = (np.argmin(np.abs(obsfwd.E()[jsmall:jlarge+1]-obsbkwd.E()[jsmall:jlarge+1]))
             +len(obsfwd.E()[:jsmall]))

        pstr = f"{k24:.2e} {gammas[j]:.2e} {jsmall} {obsfwd.R()[j]:.2e} {obsbkwd.R()[j]:.2e}"
        
        obsfwd.data[j:,:] = obsbkwd.data[j:,:]

    np.savetxt(obsfwd.observable_sname("observables",plot_format="txt"),obsfwd.data,
               fmt='\t'.join(["%13.6e"]+["%15.8e"]*5))


# Given two datasets which look at several gamma values along the coexistence line for
# an array of k24s (see image results/rasterplot-ex.pdf for what I mean), try to pin-
# point the exact coexistence line, and so trace out the differences in R,psi(R), etc
# along the coexistence line.

# To get the rastering dataset, run parameterize-coexistence.py (note that it is not
# a quick calculation, so proceed with caution to not overwrite things).


import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
from observabledata import ObservableData
from metastableregion import find_metastable

if __name__ == "__main__":



    Lambda = 27
    omega = 10
    k24s = np.linspace(0.5,1.0,num=51,endpoint=True)

    scan = {}
    scan['\\Lambda'] = str(Lambda)
    scan['\\omega']= str(omega)

    loadsuf=["K_{33}","\\Lambda","\\omega"]
    savesuf=["K_{33}","\\Lambda","\\omega"]


    obsfwd = ObservableData(["\\gamma_s","k_{24}"],scan_dir="scanforward",
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    obsbkwd = ObservableData(["\\gamma_s","k_{24}"],scan_dir="scanbackward",
                             scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    fw_coexist = np.empty([len(k24s),7],float)
    bw_coexist = np.empty([len(k24s),7],float)


    for i,k24 in enumerate(k24s):

        fw = np.where(np.isclose(obsfwd.data[:,1].T,k24),obsfwd.data.T,0).T

        fw = fw[~np.all(fw==0,axis=1)]

        bw = np.where(np.isclose(obsbkwd.data[:,1].T,k24),obsbkwd.data.T,0).T

        bw = bw[~np.all(bw==0,axis=1)]

        bw = bw[bw[:,0].argsort()]

        if len(fw) != len(bw):

            n = min(len(fw[:,0]),len(bw[:,0]))

            j = 0

            while (fw[j,0] != bw[j,0] and j < len(fw[:,0])):

                j += 1

            out_idx = np.flatnonzero(fw[j:n,0]==bw[:n-j,0])

            fw = fw[out_idx+j]

            bw = bw[out_idx]
            

        jsmall,jlarge,jeq = find_metastable(fw[:,3],bw[:,3],fw[:,2],bw[:,2])

        fw_line = []
        bw_line = []

        for k in range(7):

            fw_line.append(fw[jeq,k])
            bw_line.append(bw[jeq,k])

        fw_coexist[i,:] = fw_line
        bw_coexist[i,:] = bw_line


    np.savetxt(obsfwd.observable_sname('fwd_coexist',plot_format='txt'),fw_coexist,fmt='%13.6e')
    np.savetxt(obsbkwd.observable_sname('bkwd_coexist',plot_format='txt'),bw_coexist,fmt='%13.6e')

################################################################################
# This file can only be run if the corresponding deltazero_energy.py script    #
# has been run first! Otherwise, 'observables_Emin' type file below will not   #
# exist and the script will fail. #
################################################################################

import numpy as np
import subprocess
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams




if __name__=="__main__":



    start_time = time.time()
    
    FAILED_E = 1e300

    scan = {}
    
    loadsuf=savesuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    scan_dir = "scanforward"

    # first, load the minimum for delta = 0 case, so you know the upper bound for
    # the energy minimum.

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    run = SingleRun(rp)


    strains = np.linspace(0,0.05,num=51,endpoint=True)


    for i,u in enumerate(strains):

        if i == 0:
            
            # for the zero strain case, I need to determine what eta_eq is,
            # so I run the full 3 variable (R,eta,delta) minimization.
            
            executable = "../../../bin/full3var_onerun"
            
        else:
            
            executable = "../../../bin/delta1var_onerun"
            scan['etaguess'] = str(eta_eq/(1+u))
            scan['Rguess'] = str(R_eq/np.sqrt(1+u))

        

        # read in file name info
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
        
        # create a class to do calculations with current parameters in scan.
        run = SingleRun(rp,scan_dir=scan_dir,executable=executable)

        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file(f'observables')


        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)

        if i == 0:

            # again, if the strain is zero, then I have just determined the equilibrium
            # inverse d band spacing, which I now need to set (and do so below).

            eta_eq = etai
            R_eq = Ri

        run.concatenate_observables(None,externalparam=u)

        # now just adjust my guess for delta
        
        deltaguess = str(deltai)



        if not (np.isnan(float(deltaguess))
                or abs(float(deltaguess))<1e-5):
            scan['deltaguess'] = deltaguess
            scan['deltaupper'] = '0.818'
            
            if float(deltaguess) < 0.81:
                scan['deltalower'] = str(0.95*float(deltaguess))
            else:
                scan['deltalower'] = '0.81'

    print(f"Took {(time.time()-start_time)/3600} hours to complete.")
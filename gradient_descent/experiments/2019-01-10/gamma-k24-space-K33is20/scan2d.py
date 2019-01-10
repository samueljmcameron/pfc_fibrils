################################################################################
# This file can only be run if the corresponding deltazero_energy.py script    #
# has been run first! Otherwise, 'observables_Emin' type file below will not   #
# exist and the script will fail. #
################################################################################

import numpy as np
import subprocess
import sys
import time
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":

    start_time = time.time()
    
    FAILED_E = 1e300

    lastresortdelta = False

    k24s = np.linspace(-1,1,num=101,endpoint=True)

    k24_float = k24s[int(sys.argv[3])]
    Lambda,omega,k24 = sys.argv[1],sys.argv[2],str(k24_float)


    gammas = np.linspace(0.01,0.4,num=100,endpoint=True)

    scan = {}
    scan['k_{24}'] = k24
    scan['\\Lambda'] = Lambda
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega"]
    scan_dir = "scanforward"

    i = 0
    while (i < len(gammas)):

        gamma = gammas[i]
        
        scan['\\gamma_s']=str(gamma)

        # read in file name info
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        if i == 0:
            delta0guess = rp.params['deltaguess']
            delta0upper = rp.params['deltaupper']
            delta0lower = rp.params['deltalower']
        
        # create a class to do calculations with current parameters in scan.
        run = SingleRun(rp,scan_dir=scan_dir)


        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file('observables')

        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)


        if np.isnan(Ri) or Ri <= 0 or Ei > 0.1*FAILED_E:

            # if Ri is infinite, then the calculation failed.
            # Retry it with a different initial guess.

            # remove the current observables file, so that a new one can be written.
            run.remove_file("observables")
            if abs(float(rp.params['Rguess'])-1.0)>1e-10:
                print("Trying again with Rguess = 1.0")
                Ri = 1.0
            elif abs(float(rp.params['deltaguess'])-0)>1e-5:
                print("Last guess for Ri was 1.0, setting delta = 0.0")
                lastresortdelta = True
                deltai = 0
                Ri = 1.0
            else:
                break

        else:
            # calculation ran smoothly.
            run.concatenate_observables("\\gamma_s")
            i+= 1
        
        Rguess,etaguess,deltaguess = str(Ri),str(etai),str(deltai)

        scan['Rguess'] = Rguess
        scan['Rupper'] = str(1.5*float(Rguess))
        scan['Rlower'] = str(0.75*float(Rguess))

        if not np.isnan(float(etaguess)):
            scan['etaguess'] = etaguess
            scan['etaupper'] = str(float(etaguess)+0.03)
            scan['etalower'] = str(float(etaguess)-0.03)

        
        if (not np.isnan(float(deltaguess))
            and (abs(float(deltaguess))>1e-5
                 or lastresortdelta)):

            scan['deltaguess'] = deltaguess
            upp = np.min([float(deltaguess)+0.002,0.818])
            scan['deltaupper'] = str(upp)
            low = float(deltaguess)-0.005
            scan['deltalower'] = str(low)

        else:

            scan['deltaguess'] = delta0guess
            scan['deltaupper'] = delta0upper
            scan['deltalower'] = delta0lower

    print(f"Took {(time.time()-start_time)/3600} hours to complete.")

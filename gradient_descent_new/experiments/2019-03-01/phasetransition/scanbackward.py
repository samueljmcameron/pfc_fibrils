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

    if len(sys.argv)<4:

        user_input = input("input string of a gamma,k24,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,omega = user_input.split(',')

    else:

        gamma,k24,omega = sys.argv[1],sys.argv[2],sys.argv[3]


    Lambdas = np.linspace(0,40,num=400,endpoint=True)

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","\\omega","\\gamma_s"]
    scan_dir = "scanbackward"

    i = len(Lambdas)-1
    while (i >= 0):

        Lambda = Lambdas[i]
        
        scan['\\Lambda']=str(Lambda)

        # read in file name info
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
        
        # create a class to do calculations with current parameters in scan.
        run = SingleRun(rp,scan_dir=scan_dir)

        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file('observables')

        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)
        

        if (Ei > 0.1*FAILED_E):

            # if the energy calculation fails, this will be true.
            print('hi')
            # remove current file with observables for the current Lambda value that are higher than
            # the delta = 0 energy.
            print(Ei)
            run.remove_file("observables")

            for j,Lambda in enumerate(Lambdas[i:]):

                # write the remaining values of observables as those corresponding to the delta = 0
                # case, as non-zero d-band produces a higher energy fibril.
                scan['\\Lambda']=str(Lambda)
                rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
                run = SingleRun(rp,scan_dir=scan_dir)
                run.write_observables(E0,R0,eta0,delta0,surftwist0,"\\Lambda")

            break

        if np.isnan(Ri) or Ri <= 0:

            # if Ri is infinite, then the calculation failed.
            # Retry it with a different initial guess.

            print("Ri is NAN, trying again with Rguess = 1.0")

            # remove the current observables file, so that a new one can be written.
            run.remove_file("observables")
            if abs(float(scan['Rguess'])-1.0)>1e-10:
                Ri = 1.0
            else:
                break

        else:
            # calculation ran smoothly.
            run.concatenate_observables("\\Lambda")
            i-= 1
        
        Rguess,etaguess,deltaguess = str(Ri),str(etai),str(deltai)


        scan['Rguess'] = Rguess
        scan['Rupper'] = str(1.5*float(Rguess))
        scan['Rlower'] = str(0.75*float(Rguess))

        if not np.isnan(float(etaguess)):
            scan['etaguess'] = etaguess
            scan['etaupper'] = str(float(etaguess)+0.1)
            scan['etalower'] = str(float(etaguess)-0.02)

        if not (np.isnan(float(deltaguess))
                or abs(float(deltaguess))<1e-5):
            scan['deltaguess'] = deltaguess
            scan['deltaupper'] = '0.818'
            
            if float(deltaguess) < 0.81:
                scan['deltalower'] = str(0.95*float(deltaguess))
            else:
                scan['deltalower'] = '0.81'

    print(f"Took {(time.time()-start_time)/3600} hours to complete.")

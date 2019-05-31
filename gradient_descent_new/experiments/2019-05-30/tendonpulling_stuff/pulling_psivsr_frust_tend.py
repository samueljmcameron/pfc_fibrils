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
from psidata import PsiData
from scipy.integrate import simps

def integrand(rs,psis):

    return 2*rs*psis/(rs[-1]*rs[-1])

    

if __name__=="__main__":

    newpsiname = "frusttendon_psivsr_stress"
    newobsname = "frusttendon_observables_stress"


    start_time = time.time()

    if len(sys.argv)<5:

        user_input = input("input string of a gamma,k24,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,Lambda,omega = user_input.split(',')

    else:

        gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

    
    FAILED_E = 1e300

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\omega']= omega
    scan['\\Lambda']= Lambda
    
    loadsuf=savesuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    # first, load the minimum for delta = 0 case, so you know the upper bound for
    # the energy minimum.

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    run = SingleRun(rp)


    #strains = np.linspace(0,0.05,num=501,endpoint=True)

    strains = np.array([0.0,0.0020157,0.0034803,0.004315],float)

    #if float(Lambda) > 25.0:

    #    scan['Rguess'] = str(np.log(float(Lambda)))
    #    scan['Rupper'] = str(float(scan['Rguess'])*1.3)
    #    scan['Rlower'] = str(float(scan['Rguess'])*0.75)
    #    scan['deltaguess'] = '0.815'
    #    scan['deltalower'] = '0.81'
    #    scan['etaguess'] = '6.4'
    #    scan['etaupper'] = '6.45'

        
    for i,u in enumerate(strains):

        print(u)

        if i == 0:
            
            # for the zero strain case, I need to determine what eta_eq is,
            # so I run the full 3 variable (R,eta,delta) minimization.
            
            executable = "../../../bin/full3var_psivsr"

            strain = None
            
        else:
            
            executable = "../../../bin/delta1var_psivsr"
            scan['etaguess'] = str(eta_eq/(1+u))
            scan['Rguess'] = str(R_eq/np.sqrt(1+u))
            strain = str(u)
        

        # read in file name info
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
        
        # create a class to do calculations with current parameters in scan.
        run = SingleRun(rp,executable=executable,strain=strain)

        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file(f'observables',newname=newobsname)

        run.mv_file('psivsr',newname=newpsiname,strain=strain)

        psistuff = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,
                           name=newpsiname,strain=strain)

        rs = psistuff.r()
        psis = psistuff.psi()

        psiav = simps(integrand(rs,psis),rs)

        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,etai,deltai,surftwisti = run.get_all_observables(newobsname,str2float=True)

        if i == 0:

            # again, if the strain is zero, then I have just determined the equilibrium
            # inverse d band spacing, which I now need to set (and do so below).

            eta_eq = etai
            R_eq = Ri

        # now just adjust my guess for delta
        
        deltaguess = str(deltai)

        if np.abs(deltai)<1e-5:
            break


        if not (np.isnan(float(deltaguess))
                or abs(float(deltaguess))<1e-5):
            scan['deltaguess'] = deltaguess
            scan['deltaupper'] = '0.818'
            
            if float(deltaguess) < 0.81:
                scan['deltalower'] = str(0.95*float(deltaguess))
            else:
                scan['deltalower'] = '0.81'

    print(f"Took {(time.time()-start_time)/3600} hours to complete.")

import numpy as np
import subprocess
import sys
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":

    if len(sys.argv)<4:

        # see user_inputs.md for details on what typically goes in these inputs.
        user_input = input("input string of a k24,Lambda,omega pair, "
                           "using comma as delimiter: ")
        k24,Lambda,omega = user_input.split(',')

    else:

        k24s = np.linspace(0,1.0,num=11,endpoint=True)
        k24,Lambda,omega = str(k24s[int(sys.argv[1])]),sys.argv[2],sys.argv[3]


    gammas = np.linspace(0.02,0.2,num=19,endpoint=True)

    scan = {}
    scan['k_{24}'] = k24
    scan['\Lambda']=Lambda
    scan['\omega']=omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega"]

    for gamma in gammas:
        
        scan['\gamma_s'] = str(gamma)
        
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
        
        run = SingleRun(rp)

        run.run_exe()

        #run.mv_file('psivsr')

        run.mv_file('observables')

        
        Rguess,etaguess,deltaguess = run.get_xvals()

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


        run.concatenate_observables(["\\gamma_s"])

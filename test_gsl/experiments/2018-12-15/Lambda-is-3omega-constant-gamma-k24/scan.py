import numpy as np
import subprocess
import sys
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":

    if len(sys.argv)<3:

        user_input = input("input string of a gamma,k24 pair, "
                           "using comma as delimiter: ")
        gamma,k24 = user_input.split(',')

    else:

        gamma,k24 = sys.argv[1],sys.argv[2]


    omegas = np.linspace(19,30,num=12,endpoint=True)
    Lambdas = 3*omegas

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = str(gamma)

    
    loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","d_0","\\gamma_s"]

    for i,Lambda in enumerate(Lambdas):

        omega = omegas[i]
        
        scan['\\Lambda']=str(Lambda)
        scan['\\omega']=str(omega)
        
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


        run.concatenate_observables(["\\Lambda","\\omega"])

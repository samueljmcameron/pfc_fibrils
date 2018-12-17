import numpy as np
import subprocess
import sys
sys.path.append('../../modules_gammak24/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":

    if len(sys.argv)<4:

        user_input = input("input string of a gamma,k24,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,omega = user_input.split(',')

    else:

        gamma,k24,omega = sys.argv[1],sys.argv[2],sys.argv[3]


    Lambdas = np.linspace(0,30,num=31,endpoint=True)

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","d_0","\\omega","\\gamma_s"]


    # first, find minimum for delta = 0 case, so you know the upper bound for
    # the energy minimum.
    
    scan['\\omega']='0'
    scan['\\Lambda']='0'

    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)

    run = SingleRun(rp)

    run.run_exe()

    run.mv_file('observables',newname='observables_Emin')

    E0,R0,eta0,delta0,surftwist0 = run.get_all_observables('observables_Emin',
                                                           str2float=True)

    scan['\\omega']=omega

    for i,Lambda in enumerate(Lambdas):
        
        scan['\\Lambda']=str(Lambda)
        
        rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
        
        run = SingleRun(rp)

        run.run_exe()

        #run.mv_file('psivsr')

        run.mv_file('observables')

        Ei,Ri,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)
        
        if (abs((Ei-E0)/E0) < 1e-7):
            for j,Lambda in enumerate(Lambdas[i:]):
                scan['\\Lambda']=str(Lambda)
                rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=savesuf)
                run = SingleRun(rp)
                run.write_observables(E0,R0,eta0,delta0,surftwist0,"\\Lambda")

            break
        
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


        run.concatenate_observables("\\Lambda")

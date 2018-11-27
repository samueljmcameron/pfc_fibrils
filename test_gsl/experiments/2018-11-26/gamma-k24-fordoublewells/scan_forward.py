import numpy as np
import subprocess
import sys
sys.path.append('../../local_packages/')
from multirun import MultiRun

if __name__=="__main__":

    # see user_inputs.md for details on what typically goes in these inputs.
    #    o_L_input = input("input string of omega Lambda values, "
    #                  "using comma as delimiter: ")

    #omegas_Lambdas = np.array(o_L_input.split(','),float)

    omegas_Lambdas = np.linspace(21.5,21.5,num=1,endpoint=True)

    scan = {}

    for o_L in omegas_Lambdas:
        
        scan['\omega'] = str(o_L)
        scan['\Lambda'] = str(o_L)
        
        run = MultiRun("data/input.dat",scan)

        run.run_exe()

        run.mv_file('psivsr')

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


        run.concatenate_observables_with_direction('forward')



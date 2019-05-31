import numpy as np
import sys

sys.path.append('../../scripts/')
from doublewell import DoubleWell
from readparams import ReadParams

if __name__=="__main__":

    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    savesuf = loadsuf

    scan = {}

    scan['R0'] = str(4.29777048e-01)
    scan['R1'] = str(5.69248055e-01)
    scan['eta0'] = str(6.31908402e+00)
    scan['eta1'] = str(6.32575563e+00)
    scan['delta0'] = str(8.05969564e-01)
    scan['delta1'] = str(8.08987550e-01)
    scan['k_{24}'] = str(0.43)
    scan['\\gamma_s'] = str(8.39256250e-02)
    scan['\\Lambda'] = str(27.0)
    scan['\\omega'] = str(10.0)

    rp = ReadParams(loadsuf=loadsuf,savesuf=savesuf,scan=scan)

    dw = DoubleWell(rp);

    dw.run_exe()

    dw.mv_file("Evst")

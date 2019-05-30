import numpy as np
import sys

sys.path.append('../../scripts/')
from doublewell import DoubleWell
from readparams import ReadParams

if __name__=="__main__":

    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    savesuf = loadsuf

    rp = ReadParams(loadsuf=loadsuf,savesuf=savesuf)

    dw = DoubleWell(rp);

    dw.run_exe()

    dw.mv_file("Evst")

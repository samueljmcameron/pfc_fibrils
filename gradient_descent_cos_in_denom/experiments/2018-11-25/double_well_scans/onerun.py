import numpy as np
import sys

sys.path.append('../../local_packages/')
from doublewell import DoubleWell

if __name__=="__main__":

    dw = DoubleWell();

    dw.run_exe()

    dw.mv_file()

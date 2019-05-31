import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

def single_E_calc(gamma,scan,loadsuf,savesuf,scan_dir):

    scan['\\gamma_s'] = str(gamma)

    if scan_dir == "scanforward":

        Rguess0 = 0.15
        Rlower0 = 0.1
        Rupper0 = 0.2

        etaguess0 = 6.31
        etalower0 = 6.29
        etaupper0 = 6.33

        deltaguess0 = 0.78
        deltalower0 = 0.77
        deltaupper0 = 0.81

    else:

        Rguess0 = 2.5
        Rlower0 = 2.4
        Rupper0 = 2.7

        etaguess0 = 6.4
        etalower0 = 6.38
        etaupper0 = 6.42

        deltaguess0 = 0.815
        deltalower0 = 0.813
        deltaupper0 = 0.816


    scan['Rguess'] = str(Rguess0)
    scan['Rupper'] = str(Rupper0)
    scan['Rlower'] = str(Rlower0)

    scan['etaguess'] = str(etaguess0)
    scan['etaupper'] = str(etaupper0)
    scan['etalower'] = str(etalower0)

    scan['deltaguess'] = str(deltaguess0)
    scan['deltaupper'] = str(deltaupper0)
    scan['deltalower'] = str(deltalower0)


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

    run.concatenate_observables(["\\gamma_s"])

    return Ei,Ri


if __name__ == "__main__":

    Lambda = 27
    omega = 10
    K33 = 30
    k24 = float(sys.argv[1])



    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf=["K_{33}","k_{24}","\\Lambda","\\omega"]



    scan = {}
    scan['\\Lambda'] = str(Lambda)
    scan['\\omega']= str(omega)
    scan['K_{33}'] = str(K33)
    scan['k_{24}'] = str(k24)


    start_time = time.time()

    gamma0 = float(sys.argv[2])

    gamma2 = float(sys.argv[3])

    dg = 0.0005


    Ef0,Rf0 = single_E_calc(gamma0,scan,loadsuf,savesuf,"scanforward")
    Eb0,Rb0 = single_E_calc(gamma0,scan,loadsuf,savesuf,"scanbackward")

    while (np.abs(Rf0-Rb0)<1e-5):

        gamma0 += dg

        print(f"loop 1, gamma0 = {gamma0}")

        Ef0,Rf0 = single_E_calc(gamma0,scan,loadsuf,savesuf,"scanforward")
        Eb0,Rb0 = single_E_calc(gamma0,scan,loadsuf,savesuf,"scanbackward")

        
    print("finished loop 1")

    if np.abs(Ef0-Eb0)<1e-7:

        print("successfully found coexistence!")
        exit()

    elif Ef0>Eb0:
        
        print("lower bounding gamma is not small enough.")
        print("Coexistence point not bracketed from below.")
        print("exiting but FAILED!")

        exit()

    Ef2,Rf2 = single_E_calc(gamma2,scan,loadsuf,savesuf,"scanforward")
    Eb2,Rb2 = single_E_calc(gamma2,scan,loadsuf,savesuf,"scanbackward")



    while (np.abs(Rf2-Rb2)<1e-5):

        gamma2 -= dg

        print(f"loop 2, gamma2 = {gamma2}")

        Ef2,Rf2 = single_E_calc(gamma2,scan,loadsuf,savesuf,"scanforward")
        Eb2,Rb2 = single_E_calc(gamma2,scan,loadsuf,savesuf,"scanbackward")

    print("finished loop 2")

    if np.abs(Ef2-Eb2)<1e-7:

        print("successfully found coexistence!")
        exit()

    elif Eb2>Ef2:

        print("upper bounding gamma is not big enough.")
        print("Coexistence point not bracketed from above.")
        print("exiting but FAILED!")

        exit()
        
    Ef1 = 1
    Eb1 = 1000

    while(np.abs(Ef1-Eb1)>1e-7):

        gamma1 = 0.5*(gamma0+gamma2)

        Ef1,Rf1 = single_E_calc(gamma1,scan,loadsuf,savesuf,"scanforward")
        Eb1,Rb1 = single_E_calc(gamma1,scan,loadsuf,savesuf,"scanbackward")

        if Ef1<Eb1:

            Ef0 = Ef1
            Eb0 = Eb1

            gamma0 = gamma1

        else:

            Ef2 = Ef1
            Eb2 = Eb1

            gamma2 = gamma1


    print("successfully found coexistence!")

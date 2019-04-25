import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams


def f(x,a,b,c):

    return a*x*x+b*x+c


Lambda = 27
omega = 10
K33 = 30

data = np.loadtxt(f"data/_coexistence_{K33:1.4e}_{Lambda:1.4e}_"
                  +f"{omega:1.4e}.txt")


gammas = data[:,0]
k24s = data[:,1]


popt = np.polyfit(k24s,gammas,2)

k24s = np.linspace(0.5,1.0,num=51,endpoint=True)

loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
savesuf=["K_{33}","\\Lambda","\\omega"]
scan_dir = sys.argv[1]



scan = {}
scan['\\Lambda'] = str(Lambda)
scan['\\omega']= str(omega)

if scan_dir == "scanforward":

    Rguess0 = 0.1
    Rupper0 = 0.2
    Rlower0 = 0.05

else:

    Rguess0 = 1.0
    Rupper0 = 1.5
    Rlower0 = 0.75


start_time = time.time()

i = 0

while (i < len(k24s)):

    k24 = k24s[i]
    gamma = f(k24,*popt)

    print(gamma,k24)

    gnews = np.linspace(gamma-0.002,gamma+0.002,num=7,endpoint=False)

    if scan_dir == "scanbackward":

        gnews = gnews[::-1]

    scan['k_{24}'] = str(k24)

    scan['Rguess'] = str(Rguess0)
    scan['Rupper'] = str(Rupper0)
    scan['Rlower'] = str(Rlower0)

    for j,gnew in enumerate(gnews):

        scan['\\gamma_s'] = str(gnew)

        print(gnew,k24)

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




        if (np.isnan(Ri) or Ri <= 0) and gnew > 0.15:

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
            run.concatenate_observables(["\\gamma_s","k_{24}"])
            j+= 1

        Rguess,etaguess,deltaguess = str(Ri),str(etai),str(deltai)

        if not np.isnan(float(Rguess)):
            scan['Rguess'] = Rguess
            scan['Rupper'] = str(1.5*float(Rguess))
            scan['Rlower'] = str(0.9*float(Rguess))

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
                scan['deltaupper'] = str(1.2*float(deltaguess))
            else:
                scan['deltalower'] = '0.81'

        if scan_dir == "scanforward" and float(Rguess)>1.0:
            break
        elif scan_dir == "scanbackward" and float(Rguess) <0.5:
            break

    i += 1


print(f"Took {(time.time()-start_time)/3600} hours to complete.")


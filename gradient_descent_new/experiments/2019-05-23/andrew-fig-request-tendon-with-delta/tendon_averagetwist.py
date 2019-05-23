import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import seaborn as sns
sys.path.append('../../../../scripts/')
from fig_settings import configure_fig_settings
sys.path.append('../../scripts/')
from observabledata import ObservableData
from readparams import ReadParams

if __name__=="__main__":

    configure_fig_settings()

    width  = 3.37
    height = width*2.0

    colors = sns.color_palette()

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    observable_list = ["stress","surfacetwist","delta"]

    gamma = 0.04

    k24 = 0.5
    Lambda = 600
    omega = 20

    q = 4 # um ^{-1}
    K22 = 6 # pN

    types = ["linear","frustrated"]

    markertypes = ["--","-."]

    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i,observable in enumerate(observable_list):

        if (i > 0):
            ax[observable] = fig.add_subplot(3,1,i+1,sharex=ax[observable_list[0]])
        else:
            ax[observable] = fig.add_subplot(3,1,i+1)


    for i,type in enumerate(types):

        scan = {}
        scan['\\gamma_s']=gamma
        scan['k_{24}']=k24
        scan['\\omega']=omega
        scan['\\Lambda']=Lambda



        obsfwd = ObservableData(["strain","averagetwist"],scan_dir=f'scanforward{type}',
                                scan=scan,loadsuf=loadsuf,
                                savesuf=savesuf)
        Req = obsfwd.R()[0]

        print(Req)

        strains = obsfwd.data[:,0]
        stresses = np.gradient(obsfwd.E(),strains)*K22*q*q/1000 # stress in kPa
        stresses[0] = 0.0
        avgtwist = obsfwd.data[:,1]
        minindex = np.argmin(avgtwist)
        ysfwd = [stresses,avgtwist,obsfwd.delta()]
        if obsfwd.E()[0] > 1e299:
            print("bad calculation at Lambda = ",Lambda)



        for j,observable in enumerate(observable_list):


            if observable == 'surfacetwist':
                ylabel = r'$<\psi>$' + ' (' + r"$\si{\radian}$" + ')'
            elif observable == "stress":
                ylabel = r"$\sigma$" + ' (' + r"$\si{\kilo\pascal}$" + ')'
            elif observable == 'delta':
                ylabel = r'$\delta/\delta_0$'
                ysfwd[j] = ysfwd[j]/np.sqrt(2/3)
            elif len(observable) > 1:
                ylabel = fr'$\{observable}$'
            else:
                ylabel = fr'${observable}$'

            a = ysfwd[j][:]
            xs = strains[a>0]*100
            ys = a[a>0]
            xs = xs[:minindex]
            ys = ys[:minindex]

            ax[observable].plot(xs,ys,markertypes[i],color=colors[i],
                                label=f"{types[i]} tendon",lw=2)

            if observable == "stress":
                Y = np.gradient(ys,xs/100)
                print(f"initial young's modulus = {Y[2]/1000} MPa")
                print(f"max young's modulus = {np.max(Y)/1000} MPa")

            else:
                slope = np.gradient(ys,xs)[2]
                print(f"for {observable}, slope is {slope}")

            ax[observable].set_ylabel(ylabel,fontsize = 10)

            
    for observable in observable_list:

        if observable == "surfacetwist":
            ax[observable].legend(frameon=False)
        elif observable == "delta":
            ax[observable].set_xlabel(r"$\epsilon\times100\%$",fontsize = 10)
            ax[observable].set_xlim(left=0)

    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("psiavg_etc-vsstrain",plot_format="pdf"))

    plt.show()









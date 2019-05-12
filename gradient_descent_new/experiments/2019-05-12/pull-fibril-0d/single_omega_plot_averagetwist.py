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
    height = width*1.2


    gamma = sys.argv[1]
    k24 = sys.argv[2]
    Lambda = sys.argv[3]
    omega = sys.argv[4]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24
    scan['\\omega']=omega
    scan['\\Lambda']=Lambda

    colors = sns.color_palette()

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    #observable_list = ["E","R","delta","surfacetwist","stress"]
    observable_list = ["stress","surfacetwist"]#,"delta"]

    


    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i,observable in enumerate(observable_list):

        if (i > 0):
            ax[observable] = fig.add_subplot(2,1,i+1,sharex=ax[observable_list[0]])
        else:
            ax[observable] = fig.add_subplot(2,1,i+1)




    obsfwd = ObservableData(["strain","averagetwist"],scan_dir='scanforward',
                            scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)
    Req = obsfwd.R()[0]

    print(Req)

    q = float(input("q value = ?")) # um^{-1}
    K22 = 6 # pN

    strains = obsfwd.data[:,0]
    stresses = np.gradient(obsfwd.E(),strains)*K22*q*q/1000 # stress in kPa
    stresses[0] = 0.0
    ysfwd = [stresses,obsfwd.data[:,1],obsfwd.delta()]
    if obsfwd.E()[0] > 1e299:
        print("bad calculation at Lambda = ",Lambda)



    for i,observable in enumerate(observable_list):


        if observable == 'surfacetwist':
            ylabel = r'$<\psi(r)>_r$' + ' (' + r"$\si{\radian}$" + ')'
        elif observable == "stress":
            ylabel = r"$\sigma$" + ' (' + r"$\si{\kilo\pascal}$" + ')'
        elif observable == 'delta':
            ylabel = r'$\delta/\delta_0$'
            ysfwd[i] = ysfwd[i]/np.sqrt(2/3)
        elif len(observable) > 1:
            ylabel = fr'$\{observable}$'
        else:
            ylabel = fr'${observable}$'

        a = ysfwd[i][:]
        xs = strains[a>0]*100
        ys = a[a>0]
        ys = ys[xs<0.33]
        xs = xs[xs<0.33]

        ax[observable].plot(xs,ys,'.-',color=colors[0],
                            label="theory")

        if observable == "stress":
            Y = np.gradient(ys,xs/100)[2]
            print(f"young's modulus = {Y/1000} MPa")

        else:
            slope = np.gradient(ys,xs)[2]
            print(f"for {observable}, slope is {slope}")

        ax[observable].set_ylabel(ylabel,fontsize = 10)

    expdata = np.loadtxt("data/BellMeek2018.txt")
    strainpoints= expdata[:,1]*100
    tilt = expdata[:,2]*np.pi/180.0
    stresspoints = expdata[:,3]
            
    for observable in observable_list:

        if observable == "stress":
            ax[observable].plot(strainpoints,stresspoints,'k^',label="expt.")
            ax[observable].set_ylim(0,stresses[-1]+10)
            

        if observable == "surfacetwist":

            ax[observable].plot(strainpoints,tilt,'k^',label="expt.")
            slope = np.gradient(tilt,strainpoints)[2]
            ax[observable].set_xlabel(r"$\epsilon\times100\%$",fontsize = 10)
            ax[observable].legend(frameon=False)
            ax[observable].set_xlim(left=0)
            print(f"for experiment, slope is {slope}")

    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("psiavg_etc-vsstrain",plot_format="pdf"))

    plt.show()









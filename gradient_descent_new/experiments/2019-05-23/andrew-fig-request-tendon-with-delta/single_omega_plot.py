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
    height = width*1.5


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
    observable_list = ["stress","surfacetwist","delta"]




    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    q = 24.0 # um^{-1}

    K_22 = 6 # pN

    for i,observable in enumerate(observable_list):

        ax[observable] = fig.add_subplot(3,1,i+1)





    obsfwd = ObservableData(["strain","averagetwist"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                            savesuf=savesuf)
    Req = obsfwd.R()[0]/q*1000 # nm


    strains = obsfwd.data[:,0]
    stresses = np.gradient(obsfwd.E(),strains)*K_22*q*q/1e3 # in kPa
    stresses[0] = 0.0

    ysfwd = [stresses,obsfwd.surfacetwist(),obsfwd.delta()]
    if obsfwd.E()[0] > 1e299:
        print("bad calculation at Lambda = ",Lambda)



    for i,observable in enumerate(observable_list):


        if observable == 'surfacetwist':
            ylabel = r'$\psi(R)$' + ' (' + r'$\si{\radian}$' + ')'
        elif observable == "stress":
            ylabel = r"$\sigma$" + ' (' + r'$\si{\kilo\pascal}$' + ')'
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

        if observable == "stress":
            ax[observable].plot(xs,ys,'.-',color=colors[0],
                                label=rf"$\Lambda={Lambda}$")
            Y = np.gradient(ys,xs/100)[2]/1e3
        else:
            ax[observable].plot(strains*100,ysfwd[i][:],'.-',color=colors[0],
                                label=rf"$\Lambda={Lambda}$")
            slope = np.gradient(ys,xs)[2]
            print(f"for {observable}, slope is {slope}")
        ax[observable].set_ylabel(ylabel,fontsize = 10)
        ax[observable].set_xlabel(r"$\epsilon\times100\%$",fontsize = 10)


            
    for observable in observable_list:

        if observable == "stress":
            ax[observable].text(0.1,100,f"Youngs modulus\n={Y:3.0f}"
                                + r"$\si{\mega\pascal}$")
        if observable == "surfacetwist":
            ax[observable].text(2,0.26,rf"$R_{{eq}}=\num{{{Req:1.1e}}}$")
            expdata = np.loadtxt("data/BellMeek2018.txt")

            strainpoints= expdata[:,1]
            tilt = expdata[:,2]*np.pi/180
            ax[observable].plot(strainpoints*100,tilt,'k^')
            slope = np.gradient(tilt,strainpoints*100)[2]
            print(f"for experiment, slope is {slope}")


    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("psi_etc-vsstrain",plot_format="pdf"))

    plt.show()









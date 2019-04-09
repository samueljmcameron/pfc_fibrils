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


    k24 = '1.0'
    gamma = '0.12'
    omega = sys.argv[1]

    scan = {}
    scan['\\gamma_s']=gamma
    scan['k_{24}']=k24
    scan['\\omega']=omega

    #Larray = np.linspace(0,1000,num=101,endpoint=True)
    #Larray[0] = 1.0

    Larray = np.array([1.0,10.0,100.0,1000.0],float)


    #Larray = np.array([1.0,10.0,150.0,650.0,810.0],float)

    #Larray = np.linspace(0,30,num=4,endpoint=True)
    #Larray[0] = 1.0

    colors = sns.color_palette("hls",len(Larray))

    Lambdas = list(map(str,Larray))

    loadsuf = ["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]
    savesuf = ["K_{33}","k_{24}","\\omega","\\gamma_s"]

    #observable_list = ["E","R","delta","surfacetwist","stress"]
    observable_list = ["stress","delta","surfacetwist"]




    fig = plt.figure()
    fig.set_size_inches(width,height)

    ax = {}

    for i,observable in enumerate(observable_list):

        ax[observable] = fig.add_subplot(3,1,i+1)




    for js,Lambda in enumerate(Lambdas):

        scan['\\Lambda']=Lambda
        print(Lambda)

        obsfwd = ObservableData(["strain"],scan_dir='scanforward',scan=scan,loadsuf=loadsuf,
                                savesuf=savesuf)
        #obsfwd.sort_observables()
        #obsfwd.remove_duplicates()

        strains = obsfwd.data[:,0]
        stresses = np.gradient(obsfwd.E(),strains)
        stresses[0] = 0.0
        #ysfwd = [obsfwd.E(),obsfwd.R(),obsfwd.delta(),obsfwd.surfacetwist(),stresses]
        ysfwd = [stresses,obsfwd.delta(),obsfwd.surfacetwist()]
        if obsfwd.E()[0] > 1e299:
            print("bad calculation at Lambda = ",Lambda)
            continue



        for i,observable in enumerate(observable_list):


            if observable == 'surfacetwist':
                ylabel = r'$\psi(R)$'
            elif observable == "stress":
                ylabel = r"$\sigma$"
            elif observable == 'delta':
                ylabel = r'$\delta/\delta_0$'
                ysfwd[i] = ysfwd[i]/np.sqrt(2/3)
            elif len(observable) > 1:
                ylabel = fr'$\{observable}$'
            else:
                ylabel = fr'${observable}$'

            if observable == "stress":
                a = ysfwd[i][:]
                xs = strains[a>0]*100
                ys = a[a>0]
                ax[observable].plot(xs,ys,'.-',color=colors[js],
                                    label=rf"$\Lambda={Lambda}$")
            else:
                ax[observable].plot(strains*100,ysfwd[i][:],'.-',color=colors[js],
                                    label=rf"$\Lambda={Lambda}$")
            ax[observable].set_ylabel(ylabel,fontsize = 10)
            ax[observable].set_xlabel(r"$\epsilon\times100\%$",fontsize = 10)
        

            
    for observable in observable_list:

        ax[observable].set_xscale('log')
        if observable == "R":
            ax[observable].legend(frameon=False)
        if observable == "stress":
            ax[observable].set_yscale('log')
        if observable == "delta":
            ax[observable].legend(frameon=False)
        if observable == "surfacetwist":
            ax[observable].set_ylim(top=0.4)
            strainpoints= np.array([0.01,1.4,2.8,5.0],float)
            tilt = np.array([16,14,12,11],float)*np.pi/180
            ax[observable].plot(strainpoints,tilt,'k^')
    fig.subplots_adjust(left=0.2,right=0.8,bottom=0.1,top=0.95,hspace=0.05)
    fig.savefig(obsfwd.observable_sname("multiLambda-observable-vsstrain",plot_format="pdf"))

    plt.show()









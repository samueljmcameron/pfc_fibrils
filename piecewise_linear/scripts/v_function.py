import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


class vFunction(object):

    def __init__(self):
        
        return

    def v_approx(self,x_1,x_2,xi,zeta):

        x_2t2 = x_2*x_2

        x_2t3 = x_2*x_2*x_2

        x_2t4 = x_2*x_2*x_2*x_2

        x_2t5 = x_2*x_2*x_2*x_2*x_2

        x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2

        x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2

        x_1t2 = x_1*x_1

        x_1t3 = x_1*x_1*x_1

        x_1t4 = x_1*x_1*x_1*x_1

        x_1t5 = x_1*x_1*x_1*x_1*x_1

        x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1

        x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1


        #a0 = -2*np.sin(2*xi)*(x_2-x_1)

        #a1 = 2*np.sin(2*xi)*(x_2-x_1)-2*np.cos(2*xi)*(x_2**2-x_1**2)

        #a2 = 2*np.cos(2*xi)*(x_2**2-x_1**2)+4/3*np.sin(2*xi)*(x_2**3-x_1t3)

        #a3 = -4/3*np.sin(2*xi)*(x_2**3-x_1t3)+2/3*np.cos(2*xi)*(x_2t4-x_1t4)

        a0 = -2*np.sin(2*xi)*x_2+2*np.sin(2*xi)*x_1

        a1 = +(-2*np.cos(2*xi)*x_2t2+2*np.cos(2*xi)*x_1t2+2*np.sin(2*xi)*x_2-2*np.sin(2*xi)*x_1)

        a2 = +(4*np.sin(2*xi)*x_2t3*(1/3.0)-4*np.sin(2*xi)*x_1t3*(1/3.0)
               +2*np.cos(2*xi)*x_2t2-2*np.cos(2*xi)*x_1t2)

        a3 = +(2*np.cos(2*xi)*x_2t4*(1/3.0)-2*np.cos(2*xi)*x_1t4*(1/3.0)
               -4*np.sin(2*xi)*x_2t3*(1/3.0)+4*np.sin(2*xi)*x_1t3*(1/3.0))

        a4 = +(-4*np.sin(2*xi)*x_2t5*(1/15.0)+4*np.sin(2*xi)*x_1t5*(1/15.0)
                -2*np.cos(2*xi)*x_2t4*(1/3.0)+2*np.cos(2*xi)*x_1t4*(1/3.0))

        a5 = +(-4*np.cos(2*xi)*x_2t6*(1/45.0)+4*np.cos(2*xi)*x_1t6*(1/45.0)
                +4*np.sin(2*xi)*x_2t5*(1/15.0)-4*np.sin(2*xi)*x_1t5*(1/15.0))

        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def v_exact(self,x_1,x_2,xi,zeta):

        
        ans = np.cos(2*(zeta*x_2+xi))-np.cos(2*(zeta*x_1+xi))

        ans *= (1-zeta)/zeta

        return ans

    def v_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.v_exact(x_1,x_2,xi,zeta),
                          self.v_approx(x_1,x_2,xi,zeta))

        return result

    def dvdzeta_approx(self,x_1,x_2,xi,zeta):

        #a0 = 2*np.sin(2*xi)*(x_2-x_1)-2*np.cos(2*xi)*(x_2**2-x_1**2)

        #a1 = 2*(2*np.cos(2*xi)*(x_2**2-x_1**2)+4/3*np.sin(2*xi)*(x_2t3-x_1t3))

        #a2 = 3*(-4/3*np.sin(2*xi)*(x_2t3-x_1t3)+2/3*np.cos(2*xi)*(x_2t4-x_1t4))

        x_2t2 = x_2*x_2

        x_2t3 = x_2*x_2*x_2

        x_2t4 = x_2*x_2*x_2*x_2

        x_2t5 = x_2*x_2*x_2*x_2*x_2

        x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2

        x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2

        x_1t2 = x_1*x_1

        x_1t3 = x_1*x_1*x_1

        x_1t4 = x_1*x_1*x_1*x_1

        x_1t5 = x_1*x_1*x_1*x_1*x_1

        x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1

        x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1


        a0 = -2*np.cos(2*xi)*x_2t2+2*np.cos(2*xi)*x_1t2+2*np.sin(2*xi)*x_2-2*np.sin(2*xi)*x_1

        a1 = +(8*np.sin(2*xi)*x_2t3*(1/3.0)-8*np.sin(2*xi)*x_1t3*(1/3.0)
               +4*np.cos(2*xi)*x_2t2-4*np.cos(2*xi)*x_1t2)

        a2 = +(2*np.cos(2*xi)*x_2t4-2*np.cos(2*xi)*x_1t4
               -4*np.sin(2*xi)*x_2t3+4*np.sin(2*xi)*x_1t3)

        a3 = +(-16*np.sin(2*xi)*x_2t5*(1/15.0)+16*np.sin(2*xi)*x_1t5*(1/15.0)
                -8*np.cos(2*xi)*x_2t4*(1/3.0)+8*np.cos(2*xi)*x_1t4*(1/3.0))

        a4 = +(-4*np.cos(2*xi)*x_2t6*(1/9.0)+4*np.cos(2*xi)*x_1t6*(1/9.0)
                +4*np.sin(2*xi)*x_2t5*(1/3.0)-4*np.sin(2*xi)*x_1t5*(1/3.0))

        a5 = +(16*np.sin(2*xi)*x_2t7*(1/105.0)-16*np.sin(2*xi)*x_1t7*(1/105.0)
               +8*np.cos(2*xi)*x_2t6*(1/15.0)-8*np.cos(2*xi)*x_1t6*(1/15.0))


        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def dvdzeta_exact(self,x_1,x_2,xi,zeta):


        result = (1/zeta-1)*(2*x_1*np.sin(2*(zeta*x_1+xi))
                             -2*x_2*np.sin(2*(zeta*x_2+xi)))

        result += -1/zeta**2*(np.cos(2*(zeta*x_2+xi))-np.cos(2*(zeta*x_1+xi)))


        return result

    def dvdzeta_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dvdzeta_exact(x_1,x_2,xi,zeta),
                          self.dvdzeta_approx(x_1,x_2,xi,zeta))

        return result

    def dvdx_1_approx(self,x_1,xi,zeta):

        #a0 = 2*np.sin(2*xi)

        #a1 = -2*np.sin(2*xi)+4*np.cos(2*xi)*x_1

        #a2 = -4*np.cos(2*xi)*x_1-4*np.sin(2*xi)*x_1**2

        #a3 = 4*np.sin(2*xi)*x_1**2-8/3*np.cos(2*xi)*x_1t3

        x_2t2 = x_2*x_2

        x_2t3 = x_2*x_2*x_2

        x_2t4 = x_2*x_2*x_2*x_2

        x_2t5 = x_2*x_2*x_2*x_2*x_2

        x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2

        x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2

        x_1t2 = x_1*x_1

        x_1t3 = x_1*x_1*x_1

        x_1t4 = x_1*x_1*x_1*x_1

        x_1t5 = x_1*x_1*x_1*x_1*x_1

        x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1

        x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1


        a0 = 2*np.sin(2*xi)

        a1 = +(4*np.cos(2*xi)*x_1-2*np.sin(2*xi))

        a2 = +(-4*np.sin(2*xi)*x_1t2-4*np.cos(2*xi)*x_1)

        a3 =+(-8*np.cos(2*xi)*x_1t3*(1/3.0)+4*np.sin(2*xi)*x_1t2)

        a4 = +(4*np.sin(2*xi)*x_1t4*(1/3.0)+8*np.cos(2*xi)*x_1t3*(1/3.0))

        a5 = +(8*np.cos(2*xi)*x_1t5*(1/15.0)-4*np.sin(2*xi)*x_1t4*(1/3.0))

        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def dvdx_1_exact(self,x_1,xi,zeta):

        result = (1-zeta)*2*np.sin(2*(zeta*x_1+xi))
        
        return result

    def dvdx_1_full(self,x_1,xi,zeta):

        result = np.where(zeta!=0,self.dvdx_1_exact(x_1,xi,zeta),
                          self.dvdx_1_approx(x_1,xi,zeta))

        return result

    def dvdx_2_approx(self,x_2,xi,zeta):

        result = -self.dvdx_1_approx(x_2,xi,zeta)

        return result
    
    def dvdx_2_exact(self,x_2,xi,zeta):

        result = -self.dvdx_1_exact(x_2,xi,zeta)

        return result

    def dvdx_2_full(self,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dvdx_2_exact(x_2,xi,zeta),
                          self.dvdx_2_approx(x_2,xi,zeta))

        return result

    def dvdxi_approx(self,x_1,x_2,xi,zeta):

        #a0 = -4*np.cos(2*xi)*(x_2-x_1)

        #a1 = 4*np.cos(2*xi)*(x_2-x_1)+4*np.sin(2*xi)*(x_2**2-x_1**2)

        #a2 = -4*np.sin(2*xi)*(x_2**2-x_1**2)+8/3*np.cos(2*xi)*(x_2t3-x_1t3)

        #a3 = -8/3*np.cos(2*xi)*(x_2t3-x_1t3)-4/3*np.sin(2*xi)*(x_2t4-x_1t4)

        x_2t2 = x_2*x_2

        x_2t3 = x_2*x_2*x_2

        x_2t4 = x_2*x_2*x_2*x_2

        x_2t5 = x_2*x_2*x_2*x_2*x_2

        x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2

        x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2

        x_1t2 = x_1*x_1

        x_1t3 = x_1*x_1*x_1

        x_1t4 = x_1*x_1*x_1*x_1

        x_1t5 = x_1*x_1*x_1*x_1*x_1

        x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1

        x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1


        a0 = -4*np.cos(2*xi)*x_2+4*np.cos(2*xi)*x_1

        a1 = +(4*np.sin(2*xi)*x_2t2-4*np.sin(2*xi)*x_1t2+4*np.cos(2*xi)*x_2-4*np.cos(2*xi)*x_1)

        a2 = +(8*np.cos(2*xi)*x_2t3*(1/3.0)-8*np.cos(2*xi)*x_1t3*(1/3.0)
               -4*np.sin(2*xi)*x_2t2+4*np.sin(2*xi)*x_1t2)

        a3 = +(-4*np.sin(2*xi)*x_2t4*(1/3.0)+4*np.sin(2*xi)*x_1t4*(1/3.0)
                -8*np.cos(2*xi)*x_2t3*(1/3.0)+8*np.cos(2*xi)*x_1t3*(1/3.0))

        a4 = +(-8*np.cos(2*xi)*x_2t5*(1/15.0)+8*np.cos(2*xi)*x_1t5*(1/15.0)
                +4*np.sin(2*xi)*x_2t4*(1/3.0)-4*np.sin(2*xi)*x_1t4*(1/3.0))

        a5 = +(8*np.sin(2*xi)*x_2t6*(1/45.0)-8*np.sin(2*xi)*x_1t6*(1/45.0)
               +8*np.cos(2*xi)*x_2t5*(1/15.0)-8*np.cos(2*xi)*x_1t5*(1/15.0))


        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def dvdxi_exact(self,x_1,x_2,xi,zeta):

        result = -2*(1-zeta)/zeta*(np.sin(2*(zeta*x_2+xi))-np.sin(2*(zeta*x_1+xi)))

        return result

    def dvdxi_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dvdxi_exact(x_1,x_2,xi,zeta),
                          self.dvdxi_approx(x_1,x_2,xi,zeta))

        return result
        


if __name__ == "__main__":

    vfunc = vFunction()

    R_c = 0.0355
    R_s = 0.865
    R = 0.88
    psip_c = 3.0
    psip_s = 0.0
    

    zetas = np.linspace(-0.1,0.1,num=601,endpoint=True)

    x_1 = R_s
    x_2 = R

    xis = (psip_s-zetas)*R_s+(psip_c-psip_s)*R_c

    axs = []

    fig1 = plt.figure()

    axs.append(fig1.add_subplot(1,1,1))

    axs[0].plot(zetas,vfunc.v_approx(x_1,x_2,xis,zetas),'bo',label='approx')

    axs[0].plot(zetas,vfunc.v_exact(x_1,x_2,xis,zetas),'r.',label='exact')

    axs[0].plot(zetas,vfunc.v_full(x_1,x_2,xis,zetas),'k-',label='full')

    axs[0].set_ylabel(r'$v$')

    fig2 = plt.figure()
    
    axs.append(fig2.add_subplot(1,1,1))

    axs[1].plot(zetas,vfunc.dvdzeta_approx(x_1,x_2,xis,zetas),'bo',label='approx')

    axs[1].plot(zetas,vfunc.dvdzeta_exact(x_1,x_2,xis,zetas),'r.',label='exact')

    axs[1].plot(zetas,vfunc.dvdzeta_full(x_1,x_2,xis,zetas),'k-',label='full')

    axs[1].set_ylabel(r'$\frac{dv}{d\zeta}$')

    fig3 = plt.figure()
    
    axs.append(fig3.add_subplot(1,1,1))

    axs[2].plot(zetas,vfunc.dvdx_1_approx(x_1,xis,zetas),'bo',label='approx')

    axs[2].plot(zetas,vfunc.dvdx_1_exact(x_1,xis,zetas),'r.',label='exact')

    axs[2].plot(zetas,vfunc.dvdx_1_full(x_1,xis,zetas),'k-',label='full')

    axs[2].set_ylabel(r'$\frac{dv}{dx_1}$')

    fig4 = plt.figure()

    axs.append(fig4.add_subplot(1,1,1))

    axs[3].plot(zetas,vfunc.dvdx_2_approx(x_2,xis,zetas),'bo',label='approx')

    axs[3].plot(zetas,vfunc.dvdx_2_exact(x_2,xis,zetas),'r.',label='exact')

    axs[3].plot(zetas,vfunc.dvdx_2_full(x_2,xis,zetas),'k-',label='full')

    axs[3].set_ylabel(r'$\frac{dv}{dx_2}$')


    fig5 = plt.figure()

    axs.append(fig5.add_subplot(1,1,1))

    axs[4].plot(zetas,vfunc.dvdxi_approx(x_1,x_2,xis,zetas),'bo',label='approx')

    axs[4].plot(zetas,vfunc.dvdxi_exact(x_1,x_2,xis,zetas),'r.',label='exact')

    axs[4].plot(zetas,vfunc.dvdxi_full(x_1,x_2,xis,zetas),'k-',label='full')

    axs[4].set_ylabel(r'$\frac{dv}{d\xi}$')



    for i in range(5):

        axs[i].set_xlabel(r'$\zeta$')

        axs[i].legend(frameon=False)


    plt.show()

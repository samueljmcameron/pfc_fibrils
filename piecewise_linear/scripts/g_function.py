import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class galphaFunction(object):

    def __init__(self,alpha=None,zerotol = 1e-9):

        self.alpha = alpha
        self.zerotol = zerotol
        self.valerr = "the argument x>=pi/2 in cos(x), making integral infinite."
        return

    def integralerror(err):

        print(f"Estimate of error in integral calculation is {err}.")

        return


    def galpha_approx(self,x_1,x_2,xi,zeta):

        if self.alpha == 1:


            a0 = (x_2**2-x_1**2)/2

            a1 = 2*(x_2**3-x_1**3)*np.tan(xi)/3

            a2 = (x_2**4-x_1**4)*(3*np.tan(xi)**2+1)/4

            a3 = 4*(x_2**5-x_1**5)*(4+3*np.tan(xi)**2)*np.tan(xi)/15

        elif self.alpha == 2:
            
            a0 = (x_2**2-x_1**2)/2

            a1 = 4*(x_2**3-x_1**3)*np.tan(xi)/3

            a2 = (x_2**4-x_1**4)*(1+5*np.tan(xi)**2)/2

            a3 = (x_2**5-x_1**5)*(60*np.tan(xi)**2+28)*np.tan(xi)/15

        else:

            print(f"alpha must be either 1 or two, not {self.alpha}")


        return 1/np.cos(xi)**(2*self.alpha)*(a0 + a1*zeta + a2*zeta**2 + a3*zeta**3)


    def integrand_galpha(self,u,xi,zeta):

        if zeta*u+xi >= np.pi/2-self.zerotol:

            ValueError(self.valerr)

        return u/(np.cos(zeta*u+xi)**(2*self.alpha))


    def galpha_exact(self,x_1,x_2,xi,zeta,show_err = False):

        result = quad(lambda u: self.integrand_galpha(u,xi,zeta),x_1,x_2)

        if show_err:

            self.integralerror(result[1])

        return result[0]

    
    def galpha_full(self,x_1,x_2,xi,zeta):

        result = self.galpha_exact(x_1,x_2,xi,zeta)

        return result

    def integrand_dgalphadzeta(self,u,xi,zeta):

        if zeta*u+xi >= np.pi/2-self.zerotol:

            ValueError(self.valerr)
        
        result = 2*self.alpha*u**2*np.sin(zeta*u+xi)

        return result/np.cos(zeta*u+xi)**(2*self.alpha+1)


    def dgalphadzeta_full(self,x_1,x_2,xi,zeta,show_err=False):

        result = quad(lambda u: self.integrand_dgalaphdzeta(u,xi,zeta),x_1,x_2)

        if show_err:

            self.integralerror(result[1])

        return result[0]

    def dgalphadx_1_full(self,x_1,xi,zeta):

        return -1*self.integrand_galpha(x_1,xi,zeta)

    def dgalphadx_2_full(self,x_2,xi,zeta):

        return self.integrand_galpha(x_2,xi,zeta)

    def integrand_dgalphadxi(self,u,xi,zeta):

        if zeta*u+xi >= np.pi/2-self.zerotol:

            ValueError(self.valerr)
        
        result = 2*self.alpha*u*np.sin(zeta*u+xi)

        return result/np.cos(zeta*u+xi)**(2*self.alpha+1)

    def dgalphadxi_full(self,x_1,x_2,xi,zeta):

        result = quad(lambda u: self.integrad_dgalphadxi(self,u,xi,zeta),x_1,x_2)

        if show_err:

            self.integralerror(result[1])

        return result[0]


if __name__ == "__main__":

    g1 = galphaFunction(1)

    g2 = galphaFunction(2)


    zetas = np.linspace(-1,1,num=201,endpoint=True)

    x_1 = 0.01
    x_2 = 0.1

    xi = 0.1

    ax = [[]]


    fig1 = plt.figure()

    ax[0].append(fig1.add_subplot(1,2,1))

    ax[0].append(fig1.add_subplot(1,2,2))


    ax[0][0].plot(zetas,g1.galpha_approx(x_1,x_2,xi,zetas),'bo',label ='approx')

    ax[0][0].plot(zetas,[g1.galpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax[0][0].plot(zetas,[g1.galpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax[0][0].set_ylabel(r'$g_1$')

    ax[0][0].legend(frameon=False)

    ax[0][1].plot(zetas,g2.galpha_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax[0][1].plot(zetas,[g2.galpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax[0][1].plot(zetas,[g2.galpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')


    ax[0][1].set_ylabel(r'$g_2$')

    ax[0][1].legend(frameon=False)




    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class galphaFunction(object):

    def __init__(self,alpha=None,zerotol = 1e-9):

        self.alpha = alpha
        self.zerotol = zerotol
        return

    def falpha_approx(self,x_1,x_2,xi,zeta):

        if self.alpha == 1:


            a0 = (x_2**2-x_1**2)/2

            a1 = 2*(x_2**3-x_1**3)*np.tan(xi)/3

            a2 = (x_2**4-x_1**4)*(3*np.tan(xi)**2+1)/4

            a3 = 4*(x_2**5-x_1**5)*(4+3*np.tan(xi)**2)*np.tan(xi)/15

        elif self.alpha == 2:
            
            ## insert code here

        else:

            print(f"alpha must be either 1 or two, not {self.alpha}")


        return 1/np.cos(xi)**(2*self.alpha)*(a0 + a1*zeta + a2*zeta**2 + a3*zeta**3)

    def integrand_approx(self,u,xi,zeta):

        if xi > self.zerotol:

            result = (2*zeta/self.alpha)**(2*self.alpha)*u**(2*self.alpha-1)

        else: # this function is an approximation of sin(x+a)/x, where x is near zero.

            # If a != 0, then this corresponds to psi(r)=psi*r+c near r = 0, which is
            # not a part of our model and so should throw an exception.

            print("u value of integrand cannot have lower integration limit of zero "
                  "if the intercept is non-zero.")

            result = np.nan

        return result


    def integrand_exact(self,u,xi,zeta):

        return np.sin(2*(zeta*u+xi)/self.alpha)**(2*self.alpha)/u

    def integrand_full(self,u,xi,zeta):

        result = np.where(u>self.zerotol,self.integrand_exact(u,xi,zeta),
                          self.integrand_approx(u,xi,zeta))

        return result

    def falpha_exact(self,x_1,x_2,xi,zeta,show_err = False):

        result = quad(lambda u: self.integrand_full(u,xi,zeta),x_1,x_2)

        if show_err:

            print(f"Estimate of error in integral calculation is {result[1]}.")

        return result[0]

    
    def falpha_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.falpha_exact(x_1,x_2,xi,zeta),
                          self.falpha_approx(x_1,x_2,xi,zeta))

        return result

if __name__ == "__main__":

    f1 = falphaFunction(1)

    f2 = falphaFunction(2)


    zetas = np.linspace(-1,1,num=201,endpoint=True)

    x_1 = 0.01
    x_2 = 0.1

    xi = 0.1


    fig = plt.figure()

    ax1 = fig.add_subplot(1,2,1)

    ax2 = fig.add_subplot(1,2,2)


    ax1.plot(zetas,f1.falpha_approx(x_1,x_2,xi,zetas),'bo')

    ax1.plot(zetas,[f1.falpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.')

    ax1.plot(zetas,[f1.falpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-')

    ax2.plot(zetas,f2.falpha_approx(x_1,x_2,xi,zetas),'bo')

    ax2.plot(zetas,[f2.falpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.')

    ax2.plot(zetas,[f2.falpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-')

    plt.show()

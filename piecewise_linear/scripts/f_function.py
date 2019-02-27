import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class falphaFunction(object):

    def __init__(self,alpha=None,zerotol = 1e-9):

        self.alpha = alpha
        self.zerotol = zerotol
        return

    def falpha_approx(self,x_1,x_2,xi,zeta):

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

        if self.alpha == 1:

            if (x_1 < self.zerotol or x_2 < self.zerotol) and xi < self.zerotol:
                
                a0 = 0.0

            else:

                a0 = np.sin(2*xi)**2*np.log(x_2/x_1)

            a1 = 4*(x_2-x_1)*np.cos(2*xi)*np.sin(2*xi)

            a2 = 2*(x_2t2-x_1t2)*(np.cos(2*xi)**2-np.sin(2*xi)**2)

            a3 = -32/9.0*(x_2t3-x_1t3)*np.sin(2*xi)*np.cos(2*xi)

            a1 = 4*np.sin(2*xi)*np.cos(2*xi)*(x_2-x_1)

            a2 = +(-2*np.cos(2*xi)**2*x_1t2+2*np.cos(2*xi)**2*x_2t2+2*np.sin(2*xi)**2*x_1t2
                    -2*np.sin(2*xi)**2*x_2t2)

            a3 =-(32/9.0)*np.sin(2*xi)*np.cos(2*xi)*(-x_1t3+x_2t3)

            a4 = +(4*np.cos(2*xi)**2*x_1t4*(1/3.0)-4*np.sin(2*xi)**2*x_1t4*(1/3.0)
                   -4*np.cos(2*xi)**2*x_2t4*(1/3.0)+4*np.sin(2*xi)**2*x_2t4*(1/3.0))

            a5 = +(128/75.0)*np.sin(2*xi)*np.cos(2*xi)*(-x_1t5+x_2t5)
            

        elif self.alpha == 2:

            if (x_1 < self.zerotol or x_2 < self.zerotol) and xi < self.zerotol:
                
                a0 = 0.0

            else:
                
                a0 = np.sin(xi)**4*np.log(x_2/x_1)

            a1 = 4*(x_2-x_1)*np.sin(xi)**3*np.cos(xi)

            a2 = (x_2t2-x_1t2)*np.sin(xi)**2*(3*np.cos(xi)**2-np.sin(xi)**2)
            
            a3 = 4/3.0*(x_2t3-x_1t3)*np.sin(xi)*np.cos(xi)*(np.cos(xi)**2-5*np.sin(xi)**2)

            a1 = +4*(np.sin(xi))**3*np.cos(xi)*(x_2-x_1)

            a2 = +((np.sin(xi))**4*x_1t2-(np.sin(xi))**4*x_2t2-3*(np.sin(xi))**2*(np.cos(xi))**2*x_1t2
                   +3*(np.sin(xi))**2*(np.cos(xi))**2*x_2t2)

            a3 = +((20*(np.sin(xi))**3*np.cos(xi)*x_1t3)/9.0-(4*np.sin(xi)*(np.cos(xi))**3*x_1t3)/3.0
                   -(20*(np.sin(xi))**3*np.cos(xi)*x_2t3)/9.0+(4*np.sin(xi)*(np.cos(xi))**3*x_2t3)/3.0)

            a4 = +(-(5*(np.sin(xi))**4*x_1t4)/12.0+2*(np.sin(xi))**2*(np.cos(xi))**2*x_1t4
                    -((np.cos(xi))**4*x_1t4)/4.0+(5*(np.sin(xi))**4*x_2t4)/12.0
                    -2*(np.sin(xi))**2*(np.cos(xi))**2*x_2t4+((np.cos(xi))**4*x_2t4)/4.0)

            a5 = +(-(68*(np.sin(xi))**3*np.cos(xi)*x_1**5)/75.0+(4*np.sin(xi)*(np.cos(xi))**3*x_1**5)/5.0
                    +(68*(np.sin(xi))**3*np.cos(xi)*x_2**5)/75.0-(4*np.sin(xi)*(np.cos(xi))**3*x_2**5)/5.0)


        else:

            print(f"alpha must be either 1 or two, not {self.alpha}")


        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def integrand_f_approx(self,u,xi,zeta):

        if xi > self.zerotol:

            result = (2*zeta/self.alpha)**(2*self.alpha)*u**(2*self.alpha-1)

        else: # this function is an approximation of sin(x+a)/x, where x is near zero.

            # If a != 0, then this corresponds to psi(r)=psi*r+c near r = 0, which is
            # not a part of our model and so should throw an exception.

            ValueError("u value of integrand cannot have lower integration limit of zero "
                       "if the intercept of psi(r) is non-zero.")
            result = np.nan

        return result


    def integrand_f_exact(self,u,xi,zeta):

        return np.sin(2*(zeta*u+xi)/self.alpha)**(2*self.alpha)/u

    def integrand_f_full(self,u,xi,zeta):

        result = np.where(u>self.zerotol,self.integrand_f_exact(u,xi,zeta),
                          self.integrand_f_approx(u,xi,zeta))

        return result

    def falpha_exact(self,x_1,x_2,xi,zeta,show_err = False):

        result = quad(lambda u: self.integrand_f_full(u,xi,zeta),x_1,x_2)

        if show_err:

            print(f"Estimate of error in integral calculation is {result[1]}.")

        return result[0]

    
    def falpha_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.falpha_exact(x_1,x_2,xi,zeta),
                          self.falpha_approx(x_1,x_2,xi,zeta))

        return result

    def dfalphadzeta_approx(self,x_1,x_2,xi,zeta):

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


        if self.alpha == 1:

            a0 = 4*(x_2-x_1)*np.cos(2*xi)*np.sin(2*xi)

            a1 = 4*(x_2t2-x_1t2)*(np.cos(2*xi)**2-np.sin(2*xi)**2)

            a2 = -32/3.0*(x_2t3-x_1t3)*np.sin(2*xi)*np.cos(2*xi)


            a0 = 4*np.sin(2*xi)*np.cos(2*xi)*(x_2-x_1)

            a1 = +(-4*np.cos(2*xi)**2*x_1t2+4*np.cos(2*xi)**2*x_2t2
                    +4*np.sin(2*xi)**2*x_1t2-4*np.sin(2*xi)**2*x_2t2)

            a2 = -(32/3.0)*np.sin(2*xi)*np.cos(2*xi)*(-x_1t3+x_2t3)

            a3 = (16*np.cos(2*xi)**2*x_1t4*(1/3.0)-16*np.sin(2*xi)**2*x_1t4*(1/3.0)
                  -16*np.cos(2*xi)**2*x_2t4*(1/3.0)+16*np.sin(2*xi)**2*x_2t4*(1/3.0))
            
            a4 = (128/15.0)*np.sin(2*xi)*np.cos(2*xi)*(-x_1**5+x_2**5)
            
            a5 = (-128*np.cos(2*xi)**2*x_1**6*(1/45.0)+128*np.sin(2*xi)**2*x_1**6*(1/45.0)
                   +128*np.cos(2*xi)**2*x_2**6*(1/45.0)-128*np.sin(2*xi)**2*x_2**6*(1/45.0))
            

        elif self.alpha == 2:

            a0 = 4*(x_2-x_1)*np.sin(xi)**3*np.cos(xi)

            a1 = 2*(x_2t2-x_1t2)*np.sin(xi)**2*(3*np.cos(xi)**2-np.sin(xi)**2)
            
            a2 = 4*(x_2t3-x_1t3)*np.sin(xi)*np.cos(xi)*(np.cos(xi)**2-5*np.sin(xi)**2)

            a0 = 4*np.sin(xi)**3*np.cos(xi)*(x_2-x_1)

            a1 = +(2*np.sin(xi)**4*x_1t2-2*np.sin(xi)**4*x_2t2-6*np.sin(xi)**2*np.cos(xi)**2*x_1t2
                   +6*np.sin(xi)**2*np.cos(xi)**2*x_2t2)

            a2 = +(20*np.sin(xi)**3*np.cos(xi)*x_1t3*(1/3.0)-4*np.sin(xi)*np.cos(xi)**3*x_1t3
                   -20*np.sin(xi)**3*np.cos(xi)*x_2t3*(1/3.0)+4*np.sin(xi)*np.cos(xi)**3*x_2t3)

            a3 = +(-5*np.sin(xi)**4*x_1t4*(1/3.0)+8*np.sin(xi)**2*np.cos(xi)**2*x_1t4
                    -np.cos(xi)**4*x_1t4+5*np.sin(xi)**4*x_2t4*(1/3.0)
                    -8*np.sin(xi)**2*np.cos(xi)**2*x_2t4+np.cos(xi)**4*x_2t4)

            a4 = +(-68*np.sin(xi)**3*np.cos(xi)*x_1**5*(1/15.0)+4*np.sin(xi)*np.cos(xi)**3*x_1**5
                    +68*np.sin(xi)**3*np.cos(xi)*x_2**5*(1/15.0)-4*np.sin(xi)*np.cos(xi)**3*x_2**5)

            a5 = +(34*np.sin(xi)**4*x_1**6*(1/45.0)-64*np.sin(xi)**2*np.cos(xi)**2*x_1**6*(1/15.0)
                   +2*np.cos(xi)**4*x_1**6*(1/3.0)-34*np.sin(xi)**4*x_2**6*(1/45.0)
                   +64*np.sin(xi)**2*np.cos(xi)**2*x_2**6*(1/15.0)-2*np.cos(xi)**4*x_2**6*(1/3.0))

        else:

            print(f"alpha must be either 1 or two, not {self.alpha}")

        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5


    def dfalphadzeta_exact(self,x_1,x_2,xi,zeta):

        return 1/(zeta)*(np.sin(2/self.alpha*(zeta*x_2+xi))**(2*self.alpha)
                           -np.sin(2/self.alpha*(zeta*x_1+xi))**(2*self.alpha))

    def dfalphadzeta_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dfalphadzeta_exact(x_1,x_2,xi,zeta),
                          self.dfalphadzeta_approx(x_1,x_2,xi,zeta))

        return result


    def dfalphadx_1_full(self,x_1,xi,zeta):

        return -1*self.integrand_f_full(x_1,xi,zeta)

    def dfalphadx_2_full(self,x_2,xi,zeta):

        return self.integrand_f_full(x_2,xi,zeta)

    def dfalphadxi_approx(self,x_1,x_2,xi,zeta):

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

        if self.alpha == 1:
            
            if (x_1 < self.zerotol or x_2 < self.zerotol) and xi < self.zerotol:
                
                a0 = 0.0

            else:

                a0 = 4*np.sin(2*xi)*np.cos(2*xi)*np.log(x_2/x_1)

            a1 = 8*(x_2-x_1)*(np.cos(2*xi)**2-np.sin(2*xi)**2)

            a1 = +(8*(np.cos(2*xi))**2*(x_2-x_1)-8*(np.sin(2*xi))**2*(x_2-x_1))

            a2 = +(16*np.cos(2*xi)*x_1t2*np.sin(2*xi)-16*np.cos(2*xi)*x_2t2*np.sin(2*xi))

            a3 = +(-(64*(np.cos(2*xi))**2*(-x_1t3+x_2t3))/9.0+(64*(np.sin(2*xi))**2*(-x_1t3+x_2t3))/9.0)

            a4 = +(-(32*np.cos(2*xi)*x_1t4*np.sin(2*xi))/3.0+(32*np.cos(2*xi)*x_2t4*np.sin(2*xi))/3.0)

            a5 = +((256*(np.cos(2*xi))**2*(-x_1**5+x_2**5))/75.0-(256*(np.sin(2*xi))**2*(-x_1**5+x_2**5))/75.0)

        elif self.alpha == 2:

            if (x_1 < self.zerotol or x_2 < self.zerotol) and xi < self.zerotol:
                
                a0 = 0.0

            else:

                a0 = 4*np.sin(xi)**3*np.cos(xi)*np.log(x_2/x_1)

            a1 = 4*(x_2-x_1)*np.sin(xi)**2*(3*np.cos(xi)**2-np.sin(xi)**2)

            a1 = +(12*(np.sin(xi))**2*(np.cos(xi))**2*(x_2-x_1)-4*(np.sin(xi))**4*(x_2-x_1))

            a2 = +(10*(np.sin(xi))**3*x_1t2*np.cos(xi)-10*(np.sin(xi))**3*x_2t2*np.cos(xi)
                   -6*np.sin(xi)*(np.cos(xi))**3*x_1t2+6*np.sin(xi)*(np.cos(xi))**3*x_2t2)

            a3 = +((32*(np.sin(xi))**2*(np.cos(xi))**2*x_1t3)/3.0-(20*(np.sin(xi))**4*x_1t3)/9.0
                   -(4*(np.cos(xi))**4*x_1t3)/3.0-(32*(np.sin(xi))**2*(np.cos(xi))**2*x_2t3)/3.0
                   +(20*(np.sin(xi))**4*x_2t3)/9.0+(4*(np.cos(xi))**4*x_2t3)/3.0)

            a4 = +(-(17*(np.sin(xi))**3*x_1t4*np.cos(xi))/3.0+5*np.sin(xi)*(np.cos(xi))**3*x_1t4
                    +(17*(np.sin(xi))**3*x_2t4*np.cos(xi))/3.0-5*np.sin(xi)*(np.cos(xi))**3*x_2t4)

            a5 = +(-(128*(np.sin(xi))**2*(np.cos(xi))**2*x_1**5)/25.0+(68*(np.sin(xi))**4*x_1**5)/75.0
                    +(4*(np.cos(xi))**4*x_1**5)/5.0+(128*(np.sin(xi))**2*(np.cos(xi))**2*x_2**5)/25.0
                    -(68*(np.sin(xi))**4*x_2**5)/75.0-(4*(np.cos(xi))**4*x_2**5)/5.0)

        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3 + a4*zeta**4 + a5*zeta**5

    def integrand_dfdxi_approx(self,u,xi,zeta):

        if xi > self.zerotol:

            result = 4*(2*zeta/self.alpha)**(2*self.alpha-1)*u**(2*self.alpha-2)

        else: # this function is an approximation of sin(x+a)/x, where x is near zero.

            # If a != 0, then this corresponds to psi(r)=psi*r+c near r = 0, which is
            # not a part of our model and so should throw an exception.

            ValueError("u value of integrand cannot have lower integration limit of zero "
                       "if the intercept of psi(r) is non-zero.")

        return result

    def integrand_dfdxi_exact(self,u,xi,zeta):

        return 4*np.sin(2/self.alpha*(zeta*u+xi))**(2*self.alpha-1)*np.cos(2/self.alpha*(zeta*u+xi))/u

    def integrand_dfdxi_full(self,u,xi,zeta):

        result = np.where(u>self.zerotol,self.integrand_dfdxi_exact(u,xi,zeta),
                          self.integrand_dfdxi_approx(u,xi,zeta))

        return result

    def dfalphadxi_exact(self,x_1,x_2,xi,zeta,show_err = False):

        result = quad(lambda u: self.integrand_dfdxi_full(u,xi,zeta),x_1,x_2)

        if show_err:

            print(f"Estimate of error in integral calculation is {result[1]}.")

        return result[0]

    def dfalphadxi_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dfalphadxi_exact(x_1,x_2,xi,zeta),
                          self.dfalphadxi_approx(x_1,x_2,xi,zeta))

        return result




if __name__ == "__main__":

    f1 = falphaFunction(1)

    f2 = falphaFunction(2)


    zetas = np.linspace(-1,1,num=201,endpoint=True)

    x_1 = 0.865
    x_2 = 0.88

    xi = 0.1


    # first plot f_1 and f_2
    
    fig1 = plt.figure()

    ax1 = []
    
    ax1.append(fig1.add_subplot(1,2,1))

    ax1.append(fig1.add_subplot(1,2,2))


    ax1[0].plot(zetas,f1.falpha_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax1[0].plot(zetas,[f1.falpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
                label='exact')

    ax1[0].plot(zetas,[f1.falpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
                label='full')

    ax1[1].set_ylabel(r'$f_1$')

    ax1[1].plot(zetas,f2.falpha_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax1[1].plot(zetas,[f2.falpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
                label='exact')

    ax1[1].plot(zetas,[f2.falpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
                label='full')

    ax1[1].set_ylabel(r'$f_2$')

    ax1[1].legend(frameon=False)

    # next plot df_1/dzeta and df_2/dzeta vs zeta
    
    fig2 = plt.figure()

    ax2 = []
    
    ax2.append(fig2.add_subplot(1,2,1))

    ax2.append(fig2.add_subplot(1,2,2))


    ax2[0].plot(zetas,f1.dfalphadzeta_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax2[0].plot(zetas,[f1.dfalphadzeta_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',label='exact')

    ax2[0].plot(zetas,[f1.dfalphadzeta_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',label='full')

    ax2[0].set_ylabel(r'$\frac{\partial f_1}{\partial\zeta}$')

    ax2[1].plot(zetas,f2.dfalphadzeta_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax2[1].plot(zetas,[f2.dfalphadzeta_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',label='exact')

    ax2[1].plot(zetas,[f2.dfalphadzeta_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',label='full')

    ax2[1].set_ylabel(r'$\frac{\partial f_2}{\partial\zeta}$')
    
    ax2[1].legend(frameon=False)

    # finally, plot df_1/dxi and df_2/dxi
    
    fig3 = plt.figure()

    ax3 = []
    
    ax3.append(fig3.add_subplot(1,2,1))

    ax3.append(fig3.add_subplot(1,2,2))


    ax3[0].plot(zetas,f1.dfalphadxi_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax3[0].plot(zetas,[f1.dfalphadxi_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',label='exact')

    ax3[0].plot(zetas,[f1.dfalphadxi_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',label='full')

    ax3[1].set_ylabel(r'$\frac{\partial f_1}{\partial\xi}$')
    
    ax3[1].plot(zetas,f2.dfalphadxi_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax3[1].plot(zetas,[f2.dfalphadxi_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',label='exact')

    ax3[1].plot(zetas,[f2.dfalphadxi_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',label='full')

    ax3[1].set_ylabel(r'$\frac{\partial f_2}{\partial\xi}$')
    
    ax3[1].legend(frameon=False)
    

    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class galphaFunction(object):

    # defining the class of functions which are related to g_alpha(x_1,x_2,xi,zeta)

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


            a0 = (x_2**2-x_1**2)*np.cos(xi)**2/2

            a1 = -2*(x_2**3-x_1**3)*np.cos(xi)*np.sin(xi)/3

            a2 = -(x_2**4-x_1**4)*(np.cos(xi)**2-np.sin(xi)**2)/4

            a3 = 4*(x_2**5-x_1**5)*np.cos(xi)*np.sin(xi)/15

        elif self.alpha == 2:
            
            a0 = (x_2**2-x_1**2)*np.cos(xi)**4/2

            a1 = -4*(x_2**3-x_1**3)*np.cos(xi)**3*np.sin(xi)/3

            a2 = -(x_2**4-x_1**4)*np.cos(xi)**2*(np.cos(xi)**2-3*np.sin(xi)**2)/2

            a3 = -(x_2**5-x_1**5)*np.cos(xi)*(5*np.cos(xi)**3-12*np.sin(xi)**3)/15

        else:

            print(f"alpha must be either 1 or two, not {self.alpha}")


        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3


    def galpha_exact(self,x_1,x_2,xi,zeta):

        s1 = 2*(zeta*x_1+xi)
        s2 = 2*(zeta*x_2+xi)

        if self.alpha == 1:

            ans = 2*zeta*(x_2*np.sin(s2)-x_1*np.sin(s1))+np.cos(s2)-np.cos(s1)
            ans += 2*zeta*zeta*(x_2*x_2-x_1*x_1)

            ans /= 8*zeta*zeta

        elif self.alpha == 2:

            t1 = 4*(zeta*x_1+xi)
            t2 = 4*(zeta*x_2+xi)

            ans = 4*zeta*(x_2*np.sin(t2)-x_1*np.sin(t1))+np.cos(t2)-np.cos(t1)
            ans += 32*zeta*(x_2*np.sin(s2)-x_1*np.sin(s1))
            ans += 16*(np.cos(s2)-np.cos(s1))+24*zeta*zeta*(x_2*x_2-x_1*x_1)

            ans /= 128*zeta*zeta
            
        return ans
    
    def galpha_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.galpha_exact(x_1,x_2,xi,zeta),
                          self.galpha_approx(x_1,x_2,xi,zeta))

        return result

    def dgalphadx_1_full(self,x_1,xi,zeta):
        
        return -1*self.integrand_galpha(x_1,xi,zeta)

    def dgalphadx_2_full(self,x_2,xi,zeta):

        return self.integrand_galpha(x_2,xi,zeta)


    def dgalphadzeta_approx(self,x_1,x_2,xi,zeta):

        if self.alpha == 1:

            a0 = -2/3*(x_2**3-x_1**3)*np.cos(xi)*np.sin(xi)
            
            a1 = -1/2*(x_2**4-x_1**4)*(np.cos(xi)**2-np.sin(xi)**2)

        elif self.alpha == 2:

            a0 = -4/3*(x_2**3-x_1**3)*np.cos(xi)**3*np.sin(xi)
            
            a1 = -(x_2**4-x_1**4)*np.cos(xi)**2*(np.cos(xi)**2-np.sin(xi)**2)


        return a0 + a1*zeta
    
    def dgalphadzeta_exact(self,x_1,x_2,xi,zeta):

        p1 = zeta*x_1+xi
        p2 = zeta*x_2+xi
        
        s1 = 2*(zeta*x_1+xi)
        s2 = 2*(zeta*x_2+xi)

        q1 = -2*x_1*x_1*zeta*zeta+2*xi*xi+1
        q2 = -2*x_2*x_2*zeta*zeta+2*xi*xi+1

        if self.alpha == 1:
            
            ans = 2*zeta*(x_2*np.sin(s2)-x_1*np.sin(s1))
            ans += q2*np.cos(s2)-q1*np.cos(s1)
            ans += -4*xi*xi*(np.cos(p2)**2-np.cos(p1)**2)

            ans /= 4*zeta*zeta*zeta

            ans *= -1

        elif self.alpha == 2:

            t1 = 4*(zeta*x_1+xi)
            t2 = 4*(zeta*x_2+xi)

            r1 = -8*x_1*x_1*zeta*zeta+8*xi*xi+1
            r2 = -8*x_2*x_2*zeta*zeta+8*xi*xi+1

            ans = 4*zeta*(x_2*np.sin(t2)-x_1*np.sin(t1))+r2*np.cos(t2)-r1*np.cos(t1)
            ans += 32*zeta*(x_2*np.sin(s2)-x_1*np.sin(s1))
            ans += 16*(q2*np.cos(s2)-q1*np.cos(s1))-64*xi*xi*(np.cos(p2)**4-np.cos(p1)**4)

            ans /= 64*zeta*zeta*zeta

            ans *= -1

        return ans


    def dgalphadzeta_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dgalphadzeta_exact(x_1,x_2,xi,zeta),
                          self.dgalphadzeta_approx(x_1,x_2,xi,zeta))

        return result


    def dgalphadxi_approx(self,x_1,x_2,xi,zeta):

        if self.alpha == 1:

            a0 = -(x_2**2-x_1**2)*np.cos(xi)*np.sin(xi)

            a1 = -2/3*(x_2**3-x_1**3)*(np.cos(xi)**2-np.sin(xi)**2)

        elif self.alpha == 2:

            a0 = -2*(x_2**2-x_1**2)*np.cos(xi)**3*np.sin(xi)

            a1 = -4/3*(x_2**3-x_1**3)*np.cos(xi)**2*(np.cos(xi)**2-3*np.sin(xi)**2)

        return a0 + a1*zeta
            
    def dgalphadxi_exact(self,x_1,x_2,xi,zeta):


        p1 = zeta*x_1+xi
        p2 = zeta*x_2+xi
        
        s1 = 2*(zeta*x_1+xi)
        s2 = 2*(zeta*x_2+xi)

        if self.alpha == 1:
            
            ans = np.sin(s2)-np.sin(s1)-(s2*np.cos(s2)-s1*np.cos(s1))
            ans += 4*xi*(np.cos(p2)**2-np.cos(p1)**2)

            ans /= 4*zeta*zeta

            ans *= -1

        elif self.alpha == 2:

            t1 = 4*(zeta*x_1+xi)
            t2 = 4*(zeta*x_2+xi)

            ans = np.sin(t2)-np.sin(t1)-(t2*np.cos(t2)-t1*np.cos(t1))
            ans += 8*(np.sin(s2)-np.sin(s1))-8*(s2*np.cos(s2)-s1*np.cos(s1))
            ans += 32*xi*(np.cos(p2)**4-np.cos(p1)**4)

            ans /= 32*zeta*zeta

            ans *= -1

        return ans

    
    def dgalphadxi_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.dgalphadxi_exact(x_1,x_2,xi,zeta),
                          self.dgalphadxi_approx(x_1,x_2,xi,zeta))

        return result


    # the rest of these functions are just the original integrals, which can be
    # solved either analytically (above) or numerically (below). I'm solving them
    # numerically as a check to ensure that I did not copy anything down incorrectly
    # above.

    def integrand_galpha(self,u,xi,zeta):

        return u*np.cos(zeta*u+xi)**(2*self.alpha)


    def galpha_integral(self,x_1,x_2,xi,zeta,show_err = False):

        result = quad(lambda u: self.integrand_galpha(u,xi,zeta),x_1,x_2)

        if show_err:

            self.integralerror(result[1])

        return result[0]


    def integrand_dgalphadzeta(self,u,xi,zeta):

        result = -2*self.alpha*u**2*np.sin(zeta*u+xi)

        return result*np.cos(zeta*u+xi)**(2*self.alpha-1)


    def dgalphadzeta_integral(self,x_1,x_2,xi,zeta,show_err=False):

        result = quad(lambda u: self.integrand_dgalphadzeta(u,xi,zeta),x_1,x_2)

        if show_err:

            self.integralerror(result[1])

        return result[0]


    def integrand_dgalphadxi(self,u,xi,zeta):
        
        result = -2*self.alpha*u*np.sin(zeta*u+xi)

        return result*np.cos(zeta*u+xi)**(2*self.alpha-1)

    def dgalphadxi_integral(self,x_1,x_2,xi,zeta,show_err=False):

        result = quad(lambda u: self.integrand_dgalphadxi(u,xi,zeta),x_1,x_2)

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


    # start by just plotting g_1 and g_2 as functions of zeta
    
    ax1 = []

    
    fig1 = plt.figure()

    ax1.append(fig1.add_subplot(1,2,1))

    ax1.append(fig1.add_subplot(1,2,2))


    ax1[0].plot(zetas,g1.galpha_approx(x_1,x_2,xi,zetas),'bo',label ='approx')

    ax1[0].plot(zetas,[g1.galpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax1[0].plot(zetas,[g1.galpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax1[0].plot(zetas,[g1.galpha_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')

    ax1[0].set_ylabel(r'$g_1$')

    ax1[0].legend(frameon=False)

    ax1[1].plot(zetas,g2.galpha_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax1[1].plot(zetas,[g2.galpha_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax1[1].plot(zetas,[g2.galpha_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax1[1].plot(zetas,[g2.galpha_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')


    ax1[1].set_ylabel(r'$g_2$')

    ax1[1].legend(frameon=False)


    # now plot dg_1/dzeta and dg_2/dzeta as functions of zeta


    
    ax2 = []

    
    fig2 = plt.figure()

    ax2.append(fig2.add_subplot(1,2,1))

    ax2.append(fig2.add_subplot(1,2,2))


    ax2[0].plot(zetas,g1.dgalphadzeta_approx(x_1,x_2,xi,zetas),'bo',label ='approx')

    ax2[0].plot(zetas,[g1.dgalphadzeta_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax2[0].plot(zetas,[g1.dgalphadzeta_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax2[0].plot(zetas,[g1.dgalphadzeta_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')

    ax2[0].set_ylabel(r'$\frac{\partial g_1}{\partial\zeta}$')

    ax2[0].legend(frameon=False)

    ax2[1].plot(zetas,g2.dgalphadzeta_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax2[1].plot(zetas,[g2.dgalphadzeta_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax2[1].plot(zetas,[g2.dgalphadzeta_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax2[1].plot(zetas,[g2.dgalphadzeta_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')


    ax2[1].set_ylabel(r'$\frac{\partial g_2}{\partial\zeta}$')

    ax2[1].legend(frameon=False)

    # finally, plot dg_1/dxi and dg_2/dxi as functions of zeta

    
    ax3 = []

    
    fig3 = plt.figure()

    ax3.append(fig3.add_subplot(1,2,1))

    ax3.append(fig3.add_subplot(1,2,2))


    ax3[0].plot(zetas,g1.dgalphadxi_approx(x_1,x_2,xi,zetas),'bo',label ='approx')

    ax3[0].plot(zetas,[g1.dgalphadxi_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax3[0].plot(zetas,[g1.dgalphadxi_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax3[0].plot(zetas,[g1.dgalphadxi_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')

    ax3[0].set_ylabel(r'$\frac{\partial g_1}{\partial\xi}$')

    ax3[0].legend(frameon=False)

    ax3[1].plot(zetas,g2.dgalphadxi_approx(x_1,x_2,xi,zetas),'bo',label='approx')

    ax3[1].plot(zetas,[g2.dgalphadxi_exact(x_1,x_2,xi,zeta) for zeta in zetas],'r.',
             label='exact')

    ax3[1].plot(zetas,[g2.dgalphadxi_full(x_1,x_2,xi,zeta) for zeta in zetas],'k-',
             label='full')

    ax3[1].plot(zetas,[g2.dgalphadxi_integral(x_1,x_2,xi,zeta) for zeta in zetas],'k--',
                  label='integral')


    ax3[1].set_ylabel(r'$\frac{\partial g_2}{\partial\xi}$')

    ax3[1].legend(frameon=False)




    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


class vFunction(object):

    def __init(self):
        
        return

    def v_approx(self,x_1,x_2,xi,zeta):

        # 0th order

        a0 = -2*np.sin(2*xi)*(x_2-x_1)

        # 1st order

        a1 = 2*np.sin(2*xi)*(x_2-x_1)-2*np.cos(2*xi)*(x_2**2-x_1**2)

        # 2nd order

        a2 = 2*np.cos(2*xi)*(x_2**2-x_1**2)+4/3*np.sin(2*xi)*(x_2**3-x_1**3)

        # 3rd order

        a3 = -4/3*np.sin(2*xi)*(x_2**3-x_1**3)+2/3*np.cos(2*xi)*(x_2**4-x_1**4)


        return a0 + a1*zeta + a2*zeta**2 + a3*zeta**3

    def v_exact(self,x_1,x_2,xi,zeta):

        ans = np.cos(2*(zeta*x_2+xi))-np.cos(2*(zeta*x_1+xi))

        ans *= (1-zeta)/zeta

        return ans

    def v_full(self,x_1,x_2,xi,zeta):

        result = np.where(zeta!=0,self.v_exact(x_1,x_2,xi,zeta),
                          self.v_approx(x_1,x_2,xi,zeta))

        return result


if __name__ == "__main__":

    vfunc = vFunction()


    zetas = np.linspace(-1,1,num=201,endpoint=True)
    x_1 = 0.01
    x_2 = 0.1
    xi = 0.1

    plt.plot(zetas,vfunc.v_approx(x_1,x_2,xi,zetas),'bo')

    plt.plot(zetas,vfunc.v_exact(x_1,x_2,xi,zetas),'r.')

    plt.plot(zetas,vfunc.v_full(x_1,x_2,xi,zetas),'k-')


    plt.show()

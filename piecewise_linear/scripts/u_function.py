import numpy as np


class uFunction(object):

    def __init__(self):

        return

    
    def ufunc(self,x_1,x_2,zeta):

        return (1-zeta)*(1-zeta)*(x_2*x_2-x_1*x_1)
    

    def dudx_1(x_1,zeta):

        return -2*(1-zeta)*(1-zeta)*x_1


    def dudx_2(x_2,zeta):

        return 2*(1-zeta)*(1-zeta)*x_2


    def dudxi():

        return 0

    def dudzeta( x_1, x_2, zeta):

        return -2*zeta*(1-zeta)*(x_2*x_2-x_1*x_1)



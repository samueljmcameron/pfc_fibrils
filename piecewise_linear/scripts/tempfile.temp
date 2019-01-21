## want to re-write this into nicer code, one module per function.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def f1(a,b,psi0,psip):  # to third order in psip

    ans = -2*np.sin(2*psi0)*(b-a) # 0th order

    ans += -2*psip*np.cos(2*psi0)*(b*b-a*a) # 1st order

    ans += 4/3*psip*psip*np.sin(2*psi0)*(b*b*b-a*a*a) # second order
    
    return ans


def f2(a,b,psi0,psip): # full psip function

    print(psip.size)
    if psip.size<2 and abs(psip)<1e-10:
        ans = -np.sin(2*psi0)*(b-a)
    else:
        ans = (np.cos(2*(psip*b+psi0))-np.cos(2*(psip*a+psi0)))
        ans /= psip
    return ans


def df1(a,b,psi0,psip):  # to second order in psip

    ans = -2*np.cos(2*psi0)*(b*b-a*a)+0*psip # 0th order

    ans += 8/3*psip*np.sin(2*psi0)*(b*b*b-a*a*a) # first order

    return ans

def df2(a,b,psi0,psip): # full derivative of function
    
    ans = -1/(psip*psip)*(np.cos(2*(psip*b+psi0))-np.cos(2*(psip*a+psi0)))

    ans += -2/psip*(b*np.sin(2*(psip*b+psi0))-a*np.sin(2*(psip*a+psi0)))

    return ans


def integrand(u,psi0,psip,alpha):

    return np.sin(2/alpha*(psip*u+psi0))**(2*alpha)/u


def integral(a,b,psi0,psip,alpha):

    return quad(lambda x: integrand(x,psi0,psip,alpha),a,b)[0]

def alphais1(a,b,psi0,psip):

    x = psip/1
    y = 2*psi0/1

    if abs(a)<1e-14:
        ans = 0
    else:
        ans = np.sin(y)**2*np.log(b/a) # order zero


    ans += 4*x*(b-a)*np.cos(y)*np.sin(y)  # order one
    ans += 2*x**2*(b**2-a**2)*(np.cos(y)**2-np.sin(y)**2) # order two
    ans += -32/9*x**3*(b**3-a**3)*np.sin(y)*np.cos(y)     # order three

    return ans


def alphais2(a,b,psi0,psip):

    x = psip/2
    y = psi0

    if abs(a)<1e-14:
        ans = 0
    else:
        ans = np.sin(y)**4*np.log(b/a) # order zero


    ans += 8*x*(b-a)*np.sin(y)**3*np.cos(y) # order one
    ans += (4*x**2*(b**2-a**2)*
            (3*np.sin(y)**2*np.cos(y)**2-np.sin(y)**4)) # order two
    ans += (32/3*x**3*(b**3-a**3)*np.sin(y)*np.cos(y)*
            (np.cos(y)**2-5*np.sin(y)**2)) # order three

    return ans


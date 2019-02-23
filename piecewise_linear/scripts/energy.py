import numpy as np
from u_function import uFunction
from v_function import vFunction
from g_function import galphaFunction
from f_function import falphaFunction

def energy(R,eta,delta,R_c,R_s,psip_c,psip_s,psip_R,
           K33,Lambda,omega,k24,gamma_s):

    psi1 = (psip_c-psip_s)*R_c
    psi2 = (psip_s-psip_R)*R_s+psi1
    
    u_0 = uFunction()
    v_0 = vFunction()
    f_1 = falphaFunction(1)
    f_2 = falphaFunction(2)
    g_1 = galphaFunction(1)
    g_2 = galphaFunction(2)

    u_c = u_0.ufunc(0,R_c,psip_c)
    u_s = u_0.ufunc(R_c,R_s,psip_s)
    u_R = u_0.ufunc(R_s,R,psip_R)

    v_c = v_0.v_full(0,R_c,0,psip_c)
    v_s = v_0.v_full(R_c,R_s,psi1,psip_s)
    v_R = v_0.v_full(R_s,R,psi2,psip_R)

    f_1c = f_1.falpha_full(0,R_c,0,psip_c)
    f_1s = f_1.falpha_full(R_c,R_s,psi1,psip_s)
    f_1R = f_1.falpha_full(R_s,R,psi2,psip_R)

    f_2c = f_2.falpha_full(0,R_c,0,psip_c)
    f_2s = f_2.falpha_full(R_c,R_s,psi1,psip_s)
    f_2R = f_2.falpha_full(R_s,R,psi2,psip_R)

    g_1c = g_1.galpha_full(0,R_c,0,psip_c)
    g_1s = g_1.galpha_full(R_c,R_s,psi1,psip_s)
    g_1R = g_1.galpha_full(R_s,R,psi2,psip_R)
    
    g_2c = g_2.galpha_full(0,R_c,0,psip_c)
    g_2s = g_2.galpha_full(R_c,R_s,psi1,psip_s)
    g_2R = g_2.galpha_full(R_s,R,psi2,psip_R)

    M_PI = np.pi
    
    a1 = 0.25*(u_c+u_s+u_R);

    a1 += 0.125*(f_1c+f_1s+f_1R);

    a1 += 0.5*K33*(f_2c+f_2s+f_2R);

    a1 += 0.25*(v_c+v_s+v_R);

    a1 *= 2/(R*R);
    
    a2 = 8*M_PI*M_PI*M_PI*M_PI*R*R;


    a2 += -8*M_PI*M_PI*eta*eta*(g_1c+g_1s+g_1R);

    a2 += eta*eta*eta*eta*(g_2c+g_2s+g_2R);

    a2 *= Lambda*delta*delta/(2*R*R);


    
    a3 = 0.5*omega*delta*delta*(0.75*delta*delta-1);

    a3 += -(1+k24)*np.sin(psip_R*R+psi2)/(R*R)+2*gamma_s/R;

    piece1 = 2/(R*R)*(0.25*u_c+0.125*f_1c+0.5*K33*f_2c+0.25*v_c)
    piece1 += Lambda*delta*delta/(2*R*R)*(8*M_PI*M_PI*M_PI*M_PI*R_c*R_c
                                          -8*M_PI*M_PI*eta*eta*g_1c
                                          +eta*eta*eta*eta*g_2c)
    print("piece1 = ",piece1)

    
    piece2 = 2/(R*R)*(0.25*u_s+0.125*f_1s+0.5*K33*f_2s+0.25*v_s)
    piece2 += Lambda*delta*delta/(2*R*R)*(8*M_PI*M_PI*M_PI*M_PI*(R_s*R_s-R_c*R_c)
                                          -8*M_PI*M_PI*eta*eta*g_1s
                                          +eta*eta*eta*eta*g_2s)
    print("piece2 = ",piece2)


    piece3 = 2/(R*R)*(0.25*u_R+0.125*f_1R+0.5*K33*f_2R+0.25*v_R)
    piece3 += Lambda*delta*delta/(2*R*R)*(8*M_PI*M_PI*M_PI*M_PI*(R*R-R_s*R_s)
                                          -8*M_PI*M_PI*eta*eta*g_1R
                                          +eta*eta*eta*eta*g_2R)
    print("piece3 = ",piece3)

    
    return a1 + a2 + a3;



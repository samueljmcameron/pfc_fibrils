#include <stdio.h>
#include <stdlib.h>
#include "headerfile."


double Efunc(double *x,struct params *p)
{
  double ufunc(double x_1,double x_2,double zeta);
  double vfunc(double x_1,double x_2,double xi,double zeta);
  double f_1func(double x_1,double x_2, double xi,double zeta);
  double f_2func(double x_1,double x_2, double xi,double zeta);
  double g_1func(double x_1,double x_2,double xi,double zeta);
  double g_2func(double x_1,double x_2,double xi,double zeta);

  R = x[0];
  eta = x[1];
  delta = x[2];
  R_c = x[3];
  R_s = x[4];
  psip_c = x[5];
  psip_s = x[6];
  psip_R = x[7];


  double a1;

  a1 = 0.25*(ufunc(0,R_c,psip_c)+ufunc(R_c,R_s,psip_s)+ufunc(R_s,R,psip_R));

  a1 += 0.125*(f_1func(0,R_c,0,psip_c)+f_1func(R_c,R_s,psi1,psip_s)
	       +f_1func(R_s,R,psi2,psip_R));

  a1 += 0.5*p->K33*(f_2func(0,R_c,0,psip_c)+f_2func(R_c,R_s,psi1,psip_s)
		    +f_2func(R_s,R,psi2,psip_R));

  a1 += 0.25*(vfunc(0,R_c,0,psip_c)+vfunc(R_s,R_s,psi1,psip_s)
	      +vfunc(R_s,R,psi2,psip_R));
  
  a1 *= 2/(R*R);

  double a2;

  a2 = 16*M_PI*M_PI*M_PI*M_PI*(g_2func(0,R_c,0,psip_c)+g_2func(R_c,R_s,psi1,psip_s)
			       +g_2func(R_s,R,psi2,psip_R));

  a2 += -8*M_PI*M_PI*eta*eta*(g_1func(0,R_c,0,psip_c)+g_1func(R_c,R_s,psi1,psip_s)
			      +g_1func(R_s,R,psi2,psip_R));

  a2 += 0.5*eta*eta*eta*eta*R*R;

  a2 *= p->Lambda*delta*delta/(2*R*R);


  double a3;

  a3 = 0.5*p->omega*delta*delta*(0.75*delta*delta-1);

  a3 += -(1+p->k24)*sin(psip_R*R+psi2)/(R*R)+2*p->gamma_s/R;

  return a1 + a2 + a3;

}

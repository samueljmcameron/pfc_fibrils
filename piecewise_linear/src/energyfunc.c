#include <stdio.h>
#include <stdlib.h>
#include "headerfile.h"


double Efunc(double R,double eta, double delta,double R_c,double R_s,
	     double psip_c,double psip_s,double psip_R,struct params *p)
{
  
  double ufunc(double x_1,double x_2,double zeta);
  double vfunc(double x_1,double x_2,double xi,double zeta);
  double f_1func(double x_1,double x_2, double xi,double zeta);
  double f_2func(double x_1,double x_2, double xi,double zeta);
  double g_1func(double x_1,double x_2,double xi,double zeta);
  double g_2func(double x_1,double x_2,double xi,double zeta);

  if (R_c<0 || R_s < R_c || R < R_s) return FAILED_E;
  
  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;
  
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

  a2 = 8*M_PI*M_PI*M_PI*M_PI*R*R;

  a2 += -8*M_PI*M_PI*eta*eta*(g_1func(0,R_c,0,psip_c)+g_1func(R_c,R_s,psi1,psip_s)
			      +g_1func(R_s,R,psi2,psip_R));

  a2 += eta*eta*eta*eta*(g_2func(0,R_c,0,psip_c)+g_2func(R_c,R_s,psi1,psip_s)
			 +g_2func(R_s,R,psi2,psip_R));

  a2 *= p->Lambda*delta*delta/(2*R*R);


  double a3;

  a3 = 0.5*p->omega*delta*delta*(0.75*delta*delta-1);

  a3 += -(1+p->k24)*sin(psip_R*R+psi2)/(R*R)+2*p->gamma_s/R;

  return a1 + a2 + a3;

}

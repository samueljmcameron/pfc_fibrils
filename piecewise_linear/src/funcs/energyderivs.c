#include <stdio.h>
#include <math.h>
#include "headerfile.h"

double ufunc(double x_1,double x_2,double zeta);
double dudx_1(double x_1,double zeta);
double dudx_2(double x_2,double zeta);
double dudzeta(double x_1,double x_2,double zeta);

double vfunc(double x_1,double x_2,double xi,double zeta);
double dvdx_1(double x_1,double xi,double zeta);
double dvdx_2(double x_2,double xi,double zeta);
double dvdxi(double x_1,double x_2,double xi, double zeta);
double dvdzeta(double x_1,double x_2,double xi,double zeta);

double f_1func(double x_1,double x_2,double xi,double zeta);
double df_1dx_1(double x_1,double xi,double zeta);
double df_1dx_2(double x_2,double xi,double zeta);
double df_1dxi(double x_1,double x_2,double xi,double zeta);
double df_1dzeta(double x_1,double x_2,double xi,double zeta);

double f_2func(double x_1,double x_2,double xi,double zeta);
double df_2dx_1(double x_1,double xi,double zeta);
double df_2dx_2(double x_2,double xi,double zeta);
double df_2dxi(double x_1,double x_2,double xi,double zeta);
double df_2dzeta(double x_1,double x_2,double xi,double zeta);

double g_1func(double x_1,double x_2,double xi,double zeta);
double dg_1dx_1(double x_1,double xi,double zeta);
double dg_1dx_2(double x_2,double xi,double zeta);
double dg_1dxi(double x_1,double x_2,double xi,double zeta);
double dg_1dzeta(double x_1,double x_2,double xi,double zeta);

double g_2func(double x_1,double x_2,double xi,double zeta);
double dg_2dx_1(double x_1,double xi,double zeta);
double dg_2dx_2(double x_2,double xi,double zeta);
double dg_2dxi(double x_1,double x_2,double xi,double zeta);
double dg_2dzeta(double x_1,double x_2,double xi,double zeta);


double dEdR(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;


  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;
  
  double a1;

  a1 = 0.25*(ufunc(0,R_c,psip_c)+ufunc(R_c,R_s,psip_s)+ufunc(R_s,R,psip_R));

  a1 += 0.125*(f_1func(0,R_c,0,psip_c)+f_1func(R_c,R_s,psi1,psip_s)
	       +f_1func(R_s,R,psi2,psip_R));

  a1 += 0.5*p->K33*(f_2func(0,R_c,0,psip_c)+f_2func(R_c,R_s,psi1,psip_s)
		    +f_2func(R_s,R,psi2,psip_R));

  a1 += 0.25*(vfunc(0,R_c,0,psip_c)+vfunc(R_c,R_s,psi1,psip_s)
	      +vfunc(R_s,R,psi2,psip_R));

  a1 *= -4/(R*R*R);

  double a2;

  a2 = (0.25*dudx_2(R,psip_R)+0.125*df_1dx_2(R,psi2,psip_R)
	+0.5*p->K33*df_2dx_2(R,psi2,psip_R)+0.25*dvdx_2(R,psi2,psip_R));
  
  a2 *= 2/(R*R);

  double a3;

  a3 = 8*M_PI*M_PI*M_PI*M_PI*R*R;

  a3 += -8*M_PI*M_PI*eta*eta*(g_1func(0,R_c,0,psip_c)+g_1func(R_c,R_s,psi1,psip_s)
			      +g_1func(R_s,R,psi2,psip_R));

  a3 += eta*eta*eta*eta*(g_2func(0,R_c,0,psip_c)+g_2func(R_c,R_s,psi1,psip_s)
			 +g_2func(R_s,R,psi2,psip_R));

  a3 *= -p->Lambda*delta*delta/(R*R*R);

  double a4;

  a4 = (16*M_PI*M_PI*M_PI*M_PI*R-8*M_PI*M_PI*eta*eta*dg_1dx_2(R,psi2,psip_R)
	+eta*eta*eta*eta*dg_2dx_2(R,psi2,psip_R));

  a4 *= p->Lambda*delta*delta/(2*R*R);

  double a5;

  a5 = 2*(1+p->k24)/(R*R*R)*sin(psip_R*R+psi2);

  a5 += -psip_R*(1+p->k24)/(R*R)*cos(psip_R*R+psi2);

  a5 += -2*p->gamma_s/(R*R);

  return a1+a2+a3+a4+a5;

}


double dEdeta(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;

  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;
  
  double ans;

  ans = -16*M_PI*M_PI*eta*(g_1func(0,R_c,0,psip_c)+g_1func(R_c,R_s,psi1,psip_s)
			  +g_1func(R_s,R,psi2,psip_R));

  ans += 4*eta*eta*eta*(g_2func(0,R_c,0,psip_c)+g_2func(R_c,R_s,psi1,psip_s)
		       +g_2func(R_s,R,psi2,psip_R));

  ans *= p->Lambda*delta*delta/(2*R*R);

  return ans;

}

double dEddelta(struct params *p)
{


  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;

  
  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;
  
  double ans;

  ans = 8*M_PI*M_PI*M_PI*M_PI*R*R;

  ans += -8*M_PI*M_PI*eta*eta*(g_1func(0,R_c,0,psip_c)+g_1func(R_c,R_s,psi1,psip_s)
			       +g_1func(R_s,R,psi2,psip_R));

  ans += eta*eta*eta*eta*(g_2func(0,R_c,0,psip_c)+g_2func(R_c,R_s,psi1,psip_s)
			  +g_2func(R_s,R,psi2,psip_R));

  ans *= p->Lambda*delta/(R*R);

  ans += p->omega*delta*(1.5*delta*delta-1);

  return ans;

}

double dEdR_c(struct params *p)
{

  
  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;

  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;

  double dpsi1dR_c = psip_c-psip_s;
  double dpsi2dR_c = psip_c-psip_s;

  double a1;

  a1 = 0.25*(dudx_2(R_c,psip_c)+dudx_1(R_c,psip_s));

  a1 += 0.125*(df_1dx_2(R_c,0,psip_c)+df_1dx_1(R_c,psi1,psip_s)
	       +df_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dR_c
	       +df_1dxi(R_s,R,psi2,psip_R)*dpsi2dR_c);

  a1 += 0.5*p->K33*(df_2dx_2(R_c,0,psip_c)+df_2dx_1(R_c,psi1,psip_s)
		    +df_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dR_c
		    +df_2dxi(R_s,R,psi2,psip_R)*dpsi2dR_c);

  a1 += 0.25*(dvdx_2(R_c,0,psip_c)+dvdx_1(R_c,psi1,psip_s)
	      +dvdxi(R_c,R_s,psi1,psip_s)*dpsi1dR_c
	      +dvdxi(R_s,R,psi2,psip_R)*dpsi2dR_c);

  a1 *= 2/(R*R);

  double a2;

  a2 = -8*M_PI*M_PI*eta*eta*(dg_1dx_2(R_c,0,psip_c)+dg_1dx_1(R_c,psi1,psip_s)
			     +dg_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dR_c
			     +dg_1dxi(R_s,R,psi2,psip_R)*dpsi2dR_c);

  a2 += eta*eta*eta*eta*(dg_2dx_2(R_c,0,psip_c)+dg_2dx_1(R_c,psi1,psip_s)
			 +dg_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dR_c
			 +dg_2dxi(R_s,R,psi2,psip_R)*dpsi2dR_c);

  a2 *= p->Lambda*delta*delta/(2*R*R);

  double a3;

  a3 = -(1+p->k24)/(R*R)*cos(psip_R*R+psi2)*dpsi2dR_c;

  return a1+a2+a3;

}

double dEdR_s(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;
  
  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;

  double dpsi2dR_s = psip_s-psip_R;
  
  double a1;

  a1 = 0.25*(dudx_2(R_s,psip_s)+dudx_1(R_s,psip_R));

  a1 += 0.125*(df_1dx_2(R_s,psi1,psip_s)+df_1dx_1(R_s,psi2,psip_R)
	       +df_1dxi(R_s,R,psi2,psip_R)*dpsi2dR_s);

  a1 += 0.5*p->K33*(df_2dx_2(R_s,psi1,psip_s)+df_2dx_1(R_s,psi2,psip_R)
		    +df_2dxi(R_s,R,psi2,psip_R)*dpsi2dR_s);

  a1 += 0.25*(dvdx_2(R_s,psi1,psip_s)+dvdx_1(R_s,psi2,psip_R)
	      +dvdxi(R_s,R,psi2,psip_R)*dpsi2dR_s);

  a1 *= 2.0/(R*R);

  double a2;

  a2 = -8*M_PI*M_PI*eta*eta*(dg_1dx_2(R_s,psi1,psip_s)+dg_1dx_1(R_s,psi2,psip_R)
			     +dg_1dxi(R_s,R,psi2,psip_R)*dpsi2dR_s);

  a2 += eta*eta*eta*eta*(dg_2dx_2(R_s,psi1,psip_s)+dg_2dx_1(R_s,psi2,psip_R)
			 +dg_2dxi(R_s,R,psi2,psip_R)*dpsi2dR_s);

  a2 *= p->Lambda*delta*delta/(2*R*R);

  double a3;

  a3 = -(1+p->k24)/(R*R)*cos(psip_R*R+psi2)*dpsi2dR_s;

  return a1+a2+a3;
  
}

double dEdpsip_c(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;

  
  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;

  double dpsi1dpsip_c = R_c;
  double dpsi2dpsip_c = R_c;

  double a1;

  a1 = 0.25*dudzeta(0,R_c,psip_c);

  a1 += 0.125*(df_1dzeta(0,R_c,0,psip_c)+df_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_c
	       +df_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_c);

  a1 += 0.5*p->K33*(df_2dzeta(0,R_c,0,psip_c)+df_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_c
		    +df_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_c);

  a1 += 0.25*(dvdzeta(0,R_c,0,psip_c)+dvdxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_c
	      +dvdxi(R_s,R,psi2,psip_R)*dpsi2dpsip_c);

  a1 *= 2.0/(R*R);

  double a2;

  a2 = -8*M_PI*M_PI*eta*eta*(dg_1dzeta(0,R_c,0,psip_c)+dg_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_c
			     +dg_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_c);

  a2 += eta*eta*eta*eta*(dg_2dzeta(0,R_c,0,psip_c)+dg_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_c
			 +dg_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_c);

  a2 *= p->Lambda*delta*delta/(2*R*R);

  double a3;

  a3 = -(1+p->k24)/(R*R)*cos(psip_R*R+psi2)*dpsi2dpsip_c;

  return a1 + a2 + a3;
  
}


double dEdpsip_s(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_c = p->R_c;
  double R_s = p->R_s;
  double psip_c = p->psip_c;
  double psip_s = p->psip_s;
  double psip_R = p->psip_R;

  
  double psi1 = (psip_c-psip_s)*R_c;
  double psi2 = (psip_s-psip_R)*R_s+(psip_c-psip_s)*R_c;

  double dpsi1dpsip_s = -R_c;
  double dpsi2dpsip_s = R_s-R_c;

  double a1;

  a1 = 0.25*dudzeta(R_c,R_s,psip_s);

  a1 += 0.125*(df_1dzeta(R_c,R_s,psi1,psip_s)+df_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_s
	       +df_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_s);

  a1 += 0.5*p->K33*(df_2dzeta(R_c,R_s,psi1,psip_s)+df_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_s
		    +df_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_s);

  a1 += 0.25*(dvdzeta(R_c,R_s,psi1,psip_s)+dvdxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_s
	      +dvdxi(R_s,R,psi2,psip_R)*dpsi2dpsip_s);

  a1 *= 2.0/(R*R);

  double a2;

  a2 = -8*M_PI*M_PI*eta*eta*(dg_1dzeta(R_c,R_s,psi1,psip_s)
			     +dg_1dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_s
			     +dg_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_s);

  a2 += eta*eta*eta*eta*(dg_2dzeta(R_c,R_s,psi1,psip_s)+dg_2dxi(R_c,R_s,psi1,psip_s)*dpsi1dpsip_s
			 +dg_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_s);

  a2 *= p->Lambda*delta*delta/(2*R*R);

  double a3;

  a3 = -(1+p->k24)/(R*R)*cos(psip_R*R+psi2)*dpsi2dpsip_s;

  return a1 + a2 + a3;

}


double dEdpsip_R(struct params *p)
{

  double R = p->R;
  double eta = p->eta;
  double delta = p->delta;
  double R_s = p->R_s;
  double psip_R = p->psip_R;

  double psi2 = (p->psip_s-p->psip_R)*R_s+(p->psip_c-p->psip_s)*p->R_c;

  double dpsi2dpsip_R = -R_s;

  double a1;

  a1 = 0.25*dudzeta(R_s,R,psip_R);

  a1 += 0.125*(df_1dzeta(R_s,R,psi2,psip_R)+df_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_R);

  a1 += 0.5*p->K33*(df_2dzeta(R_s,R,psi2,psip_R)+df_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_R);

  a1 += 0.25*(dvdzeta(R_s,R,psi2,psip_R)+dvdxi(R_s,R,psi2,psip_R)*dpsi2dpsip_R);

  a1 *= 2/(R*R);

  double a2;

  a2 = -8*M_PI*M_PI/(eta*eta)*(dg_1dzeta(R_s,R,psi2,psip_R)
			       +dg_1dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_R);

  a2 += dg_2dzeta(R_s,R,psi2,psip_R)+dg_2dxi(R_s,R,psi2,psip_R)*dpsi2dpsip_R;

  a2 *= eta*eta*eta*eta*p->Lambda*delta*delta/(2*R*R);

  double a3;

  a3 = -(1+p->k24)/(R*R)*cos(psip_R*R+psi2)*(R+dpsi2dpsip_R);
  
  return a1+a2+a3;

}

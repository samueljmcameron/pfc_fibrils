#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"

void compute_rf2233b1(double K33, double Lambda,double d0,
		      double L,double eta,double delta,
		      double *r,double **y, double *rf_,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  
  rf_[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    rf_[i] = r[i]*(-(y[2][i]+0.5*sin2y/r[i])
		   +0.5*(y[2][i]+0.5*sin2y/r[i])
		   *(y[2][i]+0.5*sin2y/r[i])
		   +0.5*K33*siny*siny*siny
		   *siny/(r[i]*r[i])
		   +Lambda*delta*delta/4.0*
		   (1+sin(2*eta*L)/(2*eta*L))
		   *(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta)
		   *(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta));
  }
  return;
}

void compute_integrand1(double Lambda,double eta, double d0,
			double L, double *r, double **y,
			double *integrand1,int mpt)
{
  int i;
  double cosy;
  
  integrand1[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    cosy = cos(y[1][i]);
    integrand1[i] = r[i]*(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta);
  }
  return;
}



void compute_integrand2(double Lambda,double eta, double d0,
			double L, double *r, double **y,
			double *integrand2,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  
  integrand2[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    integrand2[i] = r[i]*((4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta)
			 *(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta));
  }
  return;
}

double E_R(double k24,double omega,double R,double L,double eta,
	   double delta,double gamma_s,double gamma_t,double *r,
	   double **y,double integration_2233b1,int mpt)
{

  double E;
  E = 2.0/(R*R)*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = E+omega*0.5*(0.75*delta*delta-1+sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
		   +delta*delta*sin(4*eta*L)/(8*eta*L))
  // add surface term tension terms
  E = E+1.0/R*(-(1+k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/R+2.0*gamma_s);  
  E = E+2*gamma_t/L;

  return E;
}

double derivEdR(double K33, double k24,double Lambda,double d0, 
		double R,double L,double eta,double delta,double gamma_s,
		double *r,double **y, double integration_2233b1,
		int mpt)
{
  double ans;
  double sin2y = sin(2*y[1][mpt]);
  double siny = sin(y[1][mpt]);
  double cosy = cos(y[1][mpt]);

  ans = -4.0/(R*R*R)*integration_2233b1;

  ans += 2.0/R*(-(y[2][mpt]+0.5*sin2y/R)+0.5*(y[2][mpt]+0.5*sin2y/R)
  		*(y[2][mpt]+0.5*sin2y/R)
  		+0.5*K33*(siny*siny*siny*siny)/(R*R));

  ans += (0.5*Lambda*delta*delta/R*(1+sin(2*eta*L)/(2*eta*L))
	  *(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta)
	  *(4*M_PI*M_PI/(d0*d0*cosy*cosy)-eta*eta));
  
  ans += 2.0/(R*R)*((1+k24)*(siny*siny/R-siny*cosy*y[2][mpt])-gamma_s);

  return ans;
}

// done up to here 2018-08-10

double derivEdeta(double R,double Lambda,double eta, 
		  double d0, double L, double *r,double **y,
		  double integration1,double integration2)
{

  double ans;

  ans = (Lambda/(2*R*R)*(cos(2*eta*L)/eta-sin(2*eta*L)/(2*eta*eta*L))
	 *integration2);
  ans += (-2*Lambda/(R*R)*(1+sin(2*eta*L)/(2*eta*L))*eta
	  *integration1);
  return ans;
}

double derivEdL(double R,double Lambda,double eta, 
		double d0, double L, double gamma_t,
		double *r,double **y,double integration2)
{
  double ans;

  ans = (Lambda/(2*R*R)*(cos(2*eta*L)/L-sin(2*eta*L)/(2*eta*L*L))
	 *integration2);
  ans -= 2*gamma_t/(L*L);

  return ans;
}

void energy_stuff(double *E, double *dEdR,double *dEdeta, double *dEdL,
		  double R, double K33, double k24, double Lambda,
		  double eta, double d0, double L, double gamma_s,
		  double gamma_t, double *r, double **y, double *rf_,
		  double *integrand1, double *integrand2,int mpt)
{
  double integration_2233b1,integration1,integration2;

  compute_rf2233b1(K33,Lambda,eta,d0,L,r,y,rf_,mpt);
  integration_2233b1 = qromb(r,rf_,mpt);

  compute_integrand1(Lambda,eta,d0,L,r,y,integrand1,mpt);
  integration1 = qromb(r,integrand1,mpt);

  compute_integrand2(Lambda,eta,d0,L,r,y,integrand2,mpt);
  integration2 = qromb(r,integrand2,mpt);


  *E = E_R(R,k24,L,gamma_s,gamma_t,r,y,integration_2233b1,mpt);

  *dEdR = derivEdR(R,K33,k24,Lambda,eta,d0,L,gamma_s,r,y,
		   integration_2233b1,mpt);

  *dEdeta = derivEdeta(R,Lambda,eta,d0,L,r,y,integration1,
		       integration2);

  *dEdL = derivEdL(R,Lambda,eta,d0,L,gamma_t,r,y,integration2);

  return;
}

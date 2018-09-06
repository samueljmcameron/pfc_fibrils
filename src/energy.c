#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"


void compute_rf2233b1(struct params *p,double *r,
		      double **y, double *rf_,int mpt);
double derivEddelta(struct params *p,double integration2);
double derivEdeta(struct params *p,double integration1,double integration2);
double derivEdR(struct params *p,double *r,double **y,
		double integration_2233b1,int mpt);
double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt);
void compute_integrand2(struct params *p,double *r,
			double **y,double *integrand2,int mpt);
void compute_integrand1(struct params *p,double *r,
			double **y,double *integrand1,int mpt);



void compute_rf2233b1(struct params *p,double *r,
		      double **y, double *rf_,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  rf_[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    rf_[i] = r[i]*(-(y[2][i]+0.5*sin2y/r[i])
		   +0.5*(y[2][i]+0.5*sin2y/r[i])
		   *(y[2][i]+0.5*sin2y/r[i])
		   +0.5*p->K33*siny*siny*siny
		   *siny/(r[i]*r[i])
		   +p->Lambda*p->delta*p->delta/4.0
		   //		   *(1+sin(2*eta*L)/(2*eta*L))
		   *(tmp/(cosy*cosy)-p->eta*p->eta)
		   *(tmp/(cosy*cosy)-p->eta*p->eta));
  }
  return;
}

void compute_integrand1(struct params *p,double *r,
			double **y,double *integrand1,int mpt)
{
  int i;
  double cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  integrand1[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    cosy = cos(y[1][i]);
    integrand1[i] = r[i]*(tmp/(cosy*cosy)-p->eta*p->eta);
  }
  return;
}

void compute_integrand2(struct params *p,double *r,
			double **y,double *integrand2,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  integrand2[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    integrand2[i] = r[i]*((tmp/(cosy*cosy)-p->eta*p->eta)
			 *(tmp/(cosy*cosy)-p->eta*p->eta));
  }
  return;
}

double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt)
{

  double E;
  E = 2.0/(p->R*p->R)*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = (E+p->Lambda*p->delta*p->delta*p->omega*0.5
       *(0.75*p->delta*p->delta-1
	 //	 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
	 //	 +delta*delta*sin(4*eta*L)/(16*eta*L)
	 ));
  // add surface term tension terms
  E = E+1.0/p->R*(-(1+p->k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/p->R+2.0*p->gamma_s);  
  //  E = E+2*gamma_t/L;

  return E;
}

double derivEdR(struct params *p,double *r,double **y,
		double integration_2233b1,int mpt)
{
  double ans;
  double sin2y = sin(2*y[1][mpt]);
  double siny = sin(y[1][mpt]);
  double cosy = cos(y[1][mpt]);
  double tmpcos = 4*M_PI*M_PI/(p->d0*p->d0*cosy*cosy);

  ans = -4.0/(p->R*p->R*p->R)*integration_2233b1;

  ans += 2.0/p->R*(-(y[2][mpt]+0.5*sin2y/p->R)
		   +0.5*(y[2][mpt]+0.5*sin2y/p->R)
		   *(y[2][mpt]+0.5*sin2y/p->R)
		   +0.5*p->K33*(siny*siny*siny*siny)/(p->R*p->R));

  ans += (0.5*p->Lambda*p->delta*p->delta/p->R
	  //	  *(1+sin(2*eta*L)/(2*eta*L))
	  *(tmpcos-p->eta*p->eta)
	  *(tmpcos-p->eta*p->eta));
  
  ans += 2.0/(p->R*p->R)*((1+p->k24)
			  *(siny*siny/p->R-siny*cosy*y[2][mpt])
			  -p->gamma_s);

  return ans;
}

double derivEdeta(struct params *p,double integration1,double integration2)
{

  double ans=0;

  //  ans = ((0.5/(R*R)*integration2+0.5*omega*(delta*delta-1))
  //	 /eta*(cos(2*eta*L)-sin(2*eta*L)/(2*eta*L)));
  ans += (-2.0/(p->R*p->R)
	  //	  *(1+sin(2*eta*L)/(2*eta*L))
	  *p->eta*integration1);

  //  ans += (omega/(8*eta)*(cos(4*eta*L)-sin(4*eta*L)/(4*eta*L)));

  ans *= p->Lambda*p->delta*p->delta;

  return ans;
}

//double derivEdL(double Lambda,double omega,double R,double eta,
//		double delta, double gamma_t,double integration2)
//{
//  double ans;
//
//  ans = ((0.5/(R*R)*integration2
//	  +0.5*omega*(delta*delta-1))
//	 *(cos(2*eta*L)-sin(2*eta*L)/(2*eta*L)));
//
//  ans += (omega/8.0*(cos(4*eta*L)-sin(4*eta*L)/(4*eta*L)));
//  ans *= Lambda*delta*delta/L;
//  ans += -2*gamma_t/(L*L);
//
//  return ans;
//}

double derivEddelta(struct params *p,double integration2)
{
  double ans;

  ans = (1.0/(p->R*p->R)
	 //	 *(1+sin(2*eta*L)/(2*eta*L))
	 *integration2);
  ans += (p->omega*p->delta*p->delta
	  *(0.75
	    //+sin(2*eta*L)/(2*eta*L)+sin(4*eta*L)/(16*eta*L)
	    ));
  ans += (p->omega*(0.75*p->delta*p->delta-1
		 //		 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
		 //		 +delta*delta*sin(4*eta*L)/(16*eta*L)
		 ));

  ans *= p->Lambda*p->delta;

  return ans;

}


bool energy_stuff(double *E, double *dEdR,double *dEdeta,
		  double *dEddelta,struct params *p,
		  double *r,double **y,double *rf_,
		  double *integrand1,double *integrand2,int mpt)
{
  bool failure = true;
  double integration_2233b1,integration1,integration2;

  compute_rf2233b1(p,r,y,rf_,mpt);
  integration_2233b1 = qromb(r,rf_,mpt,&failure);

  if (failure) {
    //    printf("failure occurred at integration_2233b1.\n");
    return false;
  }
  compute_integrand1(p,r,y,integrand1,mpt);
  integration1 = qromb(r,integrand1,mpt,&failure);

  if (failure) {
    //    printf("failure occurred at integration1.\n");
    return false;
  }

  compute_integrand2(p,r,y,integrand2,mpt);
  integration2 = qromb(r,integrand2,mpt,&failure);

  if (failure) {
    //    printf("failure occurred at integration2.\n");
    return false;
  }

  *E = E_R(p,r,y,
	   integration_2233b1,mpt);

  *dEdR = derivEdR(p,r,y,integration_2233b1,mpt);

  *dEdeta = derivEdeta(p,integration1,
		       integration2);

  *dEddelta = derivEddelta(p,integration2);

  return true;
}


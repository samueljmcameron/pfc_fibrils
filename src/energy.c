/* gradient descent derivatives for energy function */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]


// free energy density which includes only the terms which require integration //
// (i.e. the q,K22,K33, and Lambda terms) //

double f2233b1_r(struct params *p,double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp);


// calculation of integrand array to be put into qromb for integration. //

void compute_rf2233b1(struct params *p,double *x,double *r,
		      double **y, double *rf_fib,int mpt);

// calculation of energy, given that the integral has been calculated. //

double pre_integral_E(struct params *p,double *x,double *r,double **y,
		      double integration_2233b1,int mpt);

// calculation of the energy, returns true if integral was successfully
// calculated.

bool E_R(double *E,struct params *p,double *x,double *r,
	 double **y,double *rf_fib,int mpt);



double f2233b1_r(struct params *p,double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp)
{
  double ans;
  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri)+p->Lambda*x[3]*x[3]/4.0
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2])
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2]));

  return ans;
}

void compute_rf2233b1(struct params *p,double *x,double *r,
		      double **y, double *rf_fib,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  rf_fib[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    rf_fib[i] = r[i]*f2233b1_r(p,x,r[i],siny,sin2y,cosy,y[2][i],tmp);
  }
  return;
}

double pre_integral_E(struct params *p,double *x,double *r,double **y,
		      double integration_2233b1,int mpt)
{

  double E;



  E = 2.0/(x[1]*x[1])*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = (E+x[3]*x[3]*p->omega*0.5
       *(0.75*x[3]*x[3]-1
	 //	 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
	 //	 +delta*delta*sin(4*eta*L)/(16*eta*L)
	 ));
  // add surface term tension terms
  E = E+0.5+1.0/x[1]*(-(1+p->k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/x[1]+2.0*p->gamma_s);  
  //  E = E+2*gamma_t/L;

  return E;
}



bool calc_E(double *E,struct params *p,double *x,double *r,
	    double **y,double *rf_fib,int mpt)
// computes E(R) given psi(r) (i.e. y) values //
  
{
  
  bool failure = true;
  double integration_2233b1;
  double tol0 = 1e-14;
  double tol2233b1;

  // tolerances ("tol...") are to determine how large the energy or 
  // derivative terms are without the integrals (setting them to 0).
  // The reason for this is if the integral calculation performed
  // by qromb has not converged (usually because the integral error
  // is so small that round-off error becomes an issue), if the
  // size of the error term is so small relative to the rest of the
  // function that it doesn't change the functions value (up to
  // some tolerance tolblabla), then the lack of convergence can just
  // be ignored.

  
  tol2233b1 = fabs(pre_integral_E(p,r,y,0,mpt)*x[1]*x[1]/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,x,r,y,rf_fib,mpt);
  integration_2233b1 = qromb(r,rf_fib,mpt,tol2233b1,&failure);
  
  if (failure) {
    printf("failed to integrate at x = (%e,%e,%e)\n"
	   x[1],x[2],x[3],mpt);
    return false;
  }

  *E = pre_integral_E(p,x,r,y,integration_2233b1,mpt);

  return true;
}

double dEdxi(double E_p_dxi, double E_m_dxi,double dxi)
{
  return (E_p_dxi-E_m_dxi)/(2.0*dxi);
}

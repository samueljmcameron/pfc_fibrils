/* gradient descent derivatives for energy function */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]


// free energy density which includes only the terms which require integration //
// (i.e. the q,K22,K33, and Lambda terms) //

double f2233b1_r(struct params *p,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp);


// calculation of integrand array to be put into qromb for integration. //

void compute_rf2233b1(struct params *p,double *r,
		      double **y, double *rf_fib,int mpt);



// calculation of the energy (requires integral to be calculated). //

double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt);








double f2233b1_r(struct params *p,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp)
{
  double ans;
  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri)+p->Lambda*p->delta*p->delta/4.0
	 *(tmp/(cos_yi*cos_yi)-p->eta*p->eta)
	 *(tmp/(cos_yi*cos_yi)-p->eta*p->eta));

  return ans;
}

void compute_rf2233b1(struct params *p,double *r,
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
    rf_fib[i] = r[i]*f2233b1_r(p,r[i],siny,sin2y,cosy,y[2][i],tmp);
  }
  return;
}

double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt)
{

  double E;
  E = 2.0/(p->R*p->R)*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = (E+p->delta*p->delta*p->omega*0.5
       *(0.75*p->delta*p->delta-1
	 //	 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
	 //	 +delta*delta*sin(4*eta*L)/(16*eta*L)
	 ));
  // add surface term tension terms
  E = E+0.5+1.0/p->R*(-(1+p->k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/p->R+2.0*p->gamma_s);  
  //  E = E+2*gamma_t/L;

  return E;
}



void fd_dEdx(double *dEdx,double E,double E_last,double *dx,
	     int num_x)
{
  int i;
  for (i = 1; i<= num_x) {
    if (fabs(dx[i])>1e-14) dEdx[i] = (E-E_last)/dx[i];
  }
  return;
}



bool energy_prop_with_hessian(double *E, double *dEdx,struct params *p,
			      double *r,double **y,double *rf_fib,
			      double *hessian,double *dx,int mpt)
{
  bool failure = true;
  double integration_2233b1,integration1,integration2;
  double tol0 = 1e-14;
  double tol2233b1,tol1,tol2;
  double E_last = *E;
  double *dEdx_last;

  dEdx_last = malloc(sizeof(double)*num_x);
  arr_cp(dEdx_last,dEdx,num_x);

  // tolerances ("tol...") are to determine how large the energy or 
  // derivative terms are without the integrals (setting them to 0).
  // The reason for this is if the integral calculation performed
  // by qromb has not converged (usually because the integral error
  // is so small that round-off error becomes an issue), if the
  // size of the error term is so small relative to the rest of the
  // function that it doesn't change the functions value (up to
  // some tolerance tolblabla), then the lack of convergence can just
  // be ignored.
  
  tol2233b1 = fabs(E_R(p,r,y,0,mpt)*p->R*p->R/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,r,y,rf_fib,mpt);
  integration_2233b1 = qromb(r,rf_fib,mpt,tol2233b1,&failure);

  if (failure) {
    printf("tol2233b1=%e\n",tol2233b1);
    //    printf("failure occurred at integration_2233b1.\n");
    return false;
  }
  
  *E = E_R(p,r,y,integration_2233b1,mpt);

  fd_dEdx(dEdx,*E,E_last,dx,num_x);
  
  return true;
}


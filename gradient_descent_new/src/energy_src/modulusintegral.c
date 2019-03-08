/* Functions to calculate the cost function, which is the energy per unit
   volume E(x), for a given value of x. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "nrutil.h"
#include "headerfile.h"

double modulus(struct params *p)
{

  double qromb(double *x,double *y, int xlength,double tol,bool *failure);

  void compute_integrand_array(double *integrand,double *r,double **y,int n);


  
  const double tol0 = 1e-14; // tolerance for integral

  bool failure=true; // stays true if the integration fails

  double *integrand = vector(1,p->mpt);  // array to store integrand value at each point

  
  compute_integrand_array(integrand,p->r,p->y,p->mpt); // put integrand values into array

  double K = (p->Lambda*p->delta*p->delta*p->eta*p->eta*p->eta*p->eta  // compute modulus
	      *qromb(p->r,p->rf_fib,p->mpt,tol0,&failure));


  if (failure) { // if integral fails, set K to nan
    
    printf("failed to integrate the MODULUS at (R,eta,delta) = (%e,%e,%e)\n",
	   p->R,p->eta,p->delta);

    K = sqrt(-1);
  }

  
  return K;
}

void compute_integrand_array(double *integrand,double *r,double **y,int n)
{

  
  double integrand_mod(double r,double psi);

  int i;

  for (i = 1; i <= n; i++) {

    integrand[i] = integrand_mod(r[i],y[1][i]);

  };
  
  return;
  
}


double integrand_mod(double r,double psi)
{

  return r*(5*cos(psi)*cos(psi)*cos(psi)*cos(psi)-3*cos(psi)*cos(psi));

}

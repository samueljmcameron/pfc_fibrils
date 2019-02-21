/* Functions to calculate the cost function, which is the energy per unit
   volume E(x), for a given value of x. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "funcs/nrutil.h"
#include "headerfile.h"

#define LRG_NBR 1e10



double E_calc(struct params *p,double *r, double **y,double *rf_fib,int n)
/*==============================================================================

  Purpose: This function calculates E(x), by first solving the ODE for psi(r)
  using relaxation methods (with solvde_wrapper function), and then integrating
  with successful_E_count. This is all done with *mpt grid points in r and y.
  If successful_E_count fails to calculate E(x), then the whole calculation 
  (for psi(r) and E(x)) is retried  with 2*((*mpt)-1)+1 grid points. This is 
  continued until either E(x) is successfully calculated, or the maximum
  number of grid points max_mpt is reached, in which case the program exits to 
  the system.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24), as
  well as the array info (r,y,rf_fib,etc).

  ------------------------------------------------------------------------------

  Returns: Returns the value of E(x) if the calculation is successful. If the 
  calculation is unsuccessful, the form of psi(r), dpsi/dr, and r*f_fibril(r)
  are saved in a file starting with "QROMB" and then an exit status exit(1) is
  invoked, exiting to the system.

  ============================================================================*/

{

  void propagate_r(double *r, double h,int mpt);

  void buildpsi(double *r,double **y,struct params *p,int n);

  bool successful_E_count(double *E,struct params *p,double *r,
			  double **y,double *rf_fib,int n);


  

  double E;
  double h = (p->R)/(n-1);

  propagate_r(r,h,n);

  buildpsi(r,y,p,n);

  if(!successful_E_count(&E,p,r,y,rf_fib,n)) {

    //E=FAILED_E;

    // NOTE that if the calculation gets to here, then all calculations from this
    // point on will return FAILED_E, which means the derivatives will all be zero,
    // and so the minimization will exit with all derivatives being zero, but 
    // E = FAILED_E.

    printf("failed calculation!\n");

  }

  return E;
}

void buildpsi(double *r,double **y,struct params *p,int n)
{

  int i;

  int n1 = p->R_c*n/p->R;

  int n2 = p->R_s*n/p->R;
  
  double psi1 = (p->psip_c-p->psip_s)*p->R_c;
  double psi2 = (p->psip_s-p->psip_R)*p->R_s+psi1;

  for (i = 1; i<= n1; i++) {

    y[1][i] = p->psip_c*r[i];
    y[2][i] = p->psip_c;

  }

  for (i = n1+1; i <= n2; i++) {

    y[1][i] = p->psip_s*r[i]+psi1;
    y[2][i] = p->psip_s;

  }

  for (i = n2+1; i<= n; i++) {

    y[1][i] = p->psip_R*r[i]+psi2;
    y[2][i] = p->psip_R;

  }

  return;
  
}




void propagate_r(double *r, double h,int mpt)
{
  // only change r since psi, psi' are stored from last loop
  int k;
  for (k=1;k <=mpt; k++) r[k] = (k-1)*h; 
  return;
}


bool successful_E_count(double *E,struct params *p,double *r,double **y,
			double *rf_fib,int n)
/*==============================================================================

  Purpose: Given the form of psi(r) (in the array y[1..2][1..mpt]), compute the
  integrand rf_fib[1..mpt], numerically integrate this integrand, and then
  compute E(x).

  ------------------------------------------------------------------------------

  Parameters:

  *E -- pointer to variable which will store the value of E(x), if the
  calculation is successful.

  p -- This struct has all of the constant parameter info (e.g. K33, k24), as
  well as the array info (r,y,rf_fib,etc).

  ------------------------------------------------------------------------------

  Returns: Returns true if the calculation of E(x) was successful, and stores
  E(x) in *E. Returns false if E(x) integral was not successful.

  ============================================================================*/  
  
{

  void compute_rf2233b1(double *rf_fib,struct params *p,double *r,
			double **y, int n);
  double no_integral_E(struct params *p,double psiR,
		       double integration_2233b1);

  double qromb(double *x,double *y, int xlength,double tol,bool *failure);
  
  bool failure = true;
  double integration_2233b1;
  const double tol0 = 1e-14;
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
  
  tol2233b1 = fabs(no_integral_E(p,y[1][n],0)*p->R*p->R/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(rf_fib,p,r,y,n);
  integration_2233b1 = qromb(r,rf_fib,n,tol2233b1,&failure);
  
  if (failure) {
    *E = integration_2233b1;
    printf("failed to integrate at (R,eta,delta) = (%e,%e,%e)\n",
	   p->R,p->eta,p->delta);
    return false;
  }

  *E = no_integral_E(p,y[1][n],integration_2233b1);
  
  return true;
}


double no_integral_E(struct params *p,double psiR,
		     double integration_2233b1)
/*==============================================================================

  Purpose: This function computes E(x) given that the integral terms have
  already been calculated and stored in integration_2233b1.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  psiR -- This is the surface twist of the fibrils, psi(R). In this code, this
  term would be y[1][mpt].

  integration_2233b1 -- This is the value of the integral that is computed
  elsewhere.

  ------------------------------------------------------------------------------

  Returns: Returns E(x) if the integral term integration_2233b1 is the true
  integral value.

  ============================================================================*/

{

  double E;

  // first calculate bulk energy per unit length
  E = 2.0/(p->R*p->R)*integration_2233b1; 

  // add density fluctuations term
  E = (E+p->delta*p->delta*p->omega*0.5
       *(0.75*p->delta*p->delta-1));

  // add surface term tension terms
  E = E+0.5+1.0/p->R*(-(1+p->k24)*(sin(psiR)*sin(psiR))/p->R+2.0*p->gamma_s);  
  

  return E;
}


void compute_rf2233b1(double *rf_fib,struct params *p,double *r, double **y, int n)
/*==============================================================================

  Purpose: This function computes r*f_fibril(r) for all ri in r[1..n], and
  stores it into the array rf_fib[1..n].

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ============================================================================*/
{

  double f2233b1_r(struct params *p,double ri,double sin_yi,
		   double sin_2yi,double cos_yi,double yi_p);

  int i;

  double siny, sin2y,cosy;
  
  rf_fib[1] = 0;
  
  for (i = 2; i <= n; i++) {  // compute f_fibril*r

    siny = sin(y[1][i]);

    sin2y = sin(2*y[1][i]);

    cosy = cos(y[1][i]);

    rf_fib[i] = r[i]*f2233b1_r(p,r[i],siny,sin2y,cosy,
			       y[2][i]);

  }

  return;
}


double f2233b1_r(struct params *p,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p)
/*==============================================================================

  Purpose: This function computes the integrand of the energy functional (which
  is the energy density r*f_fibril(r) in the model) at the point ri.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ri -- This is the current grid point where the function is being calculated
  at.

  sin_yi, sin_2yi, cos_yi -- These are just shortcut variables for calculating
  sin(psi(ri)), sin(2*psi(ri)), and cos(psi(ri)) (used to save time instead of
  evaluating trig functions multiple times).

  yi_p -- This is the value dpsi/dr evaluated at ri.

  ------------------------------------------------------------------------------

  Returns: Returns r*f_fibril(r) at the point ri.

  ============================================================================*/
{

  double ans;

  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri)+p->Lambda*p->delta*p->delta/4.0
	 *(4*M_PI*M_PI-p->eta*p->eta*cos_yi*cos_yi)
	 *(4*M_PI*M_PI-p->eta*p->eta*cos_yi*cos_yi));

  return ans;

}

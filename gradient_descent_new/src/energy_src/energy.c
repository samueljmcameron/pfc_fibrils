/* Functions to calculate the cost function, which is the energy per unit
   volume E(x), for a given value of x. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]
#define LRG_NBR 1e10


void my_fdf(const gsl_vector *x,void *params, double *func,gsl_vector *grad)
{
  double f(const gsl_vector *x_scale,void *ps);
  void df(const gsl_vector *x,void *ps,gsl_vector *g);
  
  *func = f(x,params);
  df(x,params,grad);
  return;
}

void df(const gsl_vector *x,void *ps,gsl_vector *g)
{
  int deriv_xi(double (*f)(const gsl_vector *,void *),const gsl_vector *x,
	       int i,void *ps,double h,double *result,double *abserr);
  double f(const gsl_vector *x_scale,void *ps);
  
  double h = 1e-5;
  int i;
  double result,abserr;
  struct params *p = ps;

  for (i = 0; i < 3; i++) {
    deriv_xi(f,x,i,ps,h,&result,&abserr);
    gsl_vector_set(g,i,result);
    if (abserr >= 0.5*CONV_MIN*(1+p->Escale)) {
      if (p->Escale != 0) {
	printf("abserr is %e, but the convergence criterion "
	       "CONV_min*(1+p->Escale) is %e!\n",
	       abserr,CONV_MIN*(1+p->Escale));
      }
    }
  }
  return;
  
}


double f(const gsl_vector *x_scale,void *ps)
{
  
  double E_calc(const double *x,struct params *ps);


  
  double *x;
  double E;
  struct params *p = ps;
  
  x = vector(1,3);

  scale_backward(x_scale,x,p);
    
  E = E_calc(x,p);


  
  free_vector(x,1,3);
  return E;
}

double E_calc(const double *x,struct params *p)
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

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Returns the value of E(x) if the calculation is successful. If the 
  calculation is unsuccessful, the form of psi(r), dpsi/dr, and r*f_fibril(r)
  are saved in a file starting with "QROMB" and then an exit status exit(1) is
  invoked, exiting to the system.

  ============================================================================*/

{

  void propagate_r(double *r, double h,int mpt);

  void interpolate_array(double *r,double **y,double *r_cp,double **y_cp,int mpt);

  void copy_2_arrays(double *r_cp,double **y_cp,double *r,double **y,
		     int last_mpt);

  void solvde_wrapper(double scalv[],struct params *p,const double *x,double h,
		      bool ignore_first_y);
  
  bool successful_E_count(double *E,struct params *p,const double *x);

  void write_QROMBfailure(struct params *p,const double *x);


  
  int last_mpt = p->mpt;


  double h;
  double E;  
  double scalv[2+1];



  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values


  if (x[1] <= 0 || x[2] <= 0 || x[2] >= 8.0
      || fabs(x[3]) >= 1.0) {
    printf("x[1] = %e, x[2] = %e, x[3] = %e\n",x[1],x[2],x[3]);
    printf("something is too big or less than zero, so returning failed calculation.\n");
    return FAILED_E;
  }
  
  
  while (p->mpt <= MAX_M) {

    h = x[1]/(p->mpt-1);    // compute stepsize in r[1..mpt] 

    if (p->mpt != last_mpt) {

      printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	     x[1],x[2],x[3]);

      copy_2_arrays(p->r_cp,p->y_cp,p->r,p->y,last_mpt); // copy arrays r and y into r_cp and y_cp
      propagate_r(p->r,h,p->mpt);
      interpolate_array(p->r,p->y,p->r_cp,p->y_cp,p->mpt); // interpolate old y values (now stored in y_cp)


      last_mpt = p->mpt;
    }

    else {
      propagate_r(p->r,h,p->mpt);
      copy_2_arrays(p->r_cp,p->y_cp,p->r,p->y,p->mpt); // copy arrays r and y into r_cp and y_cp
    }

    
    solvde_wrapper(scalv,p,x,h,false);

    if(successful_E_count(&E,p,x)) return E;
    else if (E > LRG_NBR) {
      solvde_wrapper(scalv,p,x,h,true);
      if(successful_E_count(&E,p,x)) return E;
    }
    p->mpt = (p->mpt-1)*2+1;
    
  }

  // if it makes it this far, we did not successfully compute E(R)

  //  write_QROMBfailure(p,x); // save psi(r), rf_fib(r), and exit

  // NOTE that if the calculation gets to here, then all calculations from this
  // point on will return FAILED_E, which means the derivatives will all be zero,
  // and so the minimization will exit with all derivatives being zero, but 
  // E = FAILED_E.

  printf("failed calculation!\n");

  return FAILED_E;
}

void interpolate_array(double *r,double **y,double *r_cp,double **y_cp,int mpt)
/*==============================================================================
  ============================================================================*/
{

  void quick_interp(double *xcp,double **ycp,double x,double **y,int i);

  int i;
  double dy;

  for (i = 1; i <=mpt; i++) quick_interp(r_cp,y_cp,r[i],y,i);
  printf("done interpolating.\n");

  return;
}

void quick_interp(double *xcp,double **ycp,double x,double **y,int i)
// Given two points of ycp[1][:] (y11cp,x1cp) and y(12cp,x2cp), //
// and two points of ycp[2][:] (y21cp,x1cp) and (y22cp,x2cp),   //
// interpolate to determine the value of y between the two

{
  if (i%2 != 0) {
    y[1][i] = ycp[1][(i-1)/2+1];
    y[2][i] = ycp[2][(i-1)/2+1];
  }
  else {
    double slope_1, slope_2;
    slope_1 = (ycp[1][i/2+1]-ycp[1][i/2])/(xcp[i/2+1]-xcp[i/2]);
    slope_2 = (ycp[2][i/2+1]-ycp[2][i/2])/(xcp[i/2+1]-xcp[i/2]);

    y[1][i] = slope_1*(x-xcp[i/2])+ycp[1][i/2];
    y[2][i] = slope_2*(x-xcp[i/2])+ycp[2][i/2];
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

void write_QROMBfailure(struct params *p,const double *x)
/*==============================================================================
  
  Purpose: This function saves the r, psi(r), and r*f_fibril(r) values that 
  were used in attempting to calculate E(x). This function is only called if
  the E(x) calculation fails.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Exits with exit status exit(1);

  ============================================================================*/

{

  void make_f_err(char *f_err,char *err_type,int f_err_size,struct params *p,
		  const double *x);
  
  int i;
  FILE *broken;
  int f_err_size = 200;
  char f_err[f_err_size];

  make_f_err(f_err,"QROMB",f_err_size,p,x);

  printf("failed to integrate with qromb at x = (%e,%e,%e).\n",x[1],x[2],x[3]);
  printf("saving psi(r) shape, rf_fib(r), and exiting to system.\n");

  broken = fopen(f_err,"w");


  for (i = 1; i<=p->mpt; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",p->r[i],p->y[1][i],p->y[2][i],
	    p->rf_fib[i]);
  }
  fclose(broken);
  exit(1);
  return;
}




bool successful_E_count(double *E,struct params *p,const double *x)
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

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Returns true if the calculation of E(x) was successful, and stores
  E(x) in *E. Returns false if E(x) integral was not successful.

  ============================================================================*/  
  
{

  void compute_rf2233b1(struct params *p,const double *x);
  double no_integral_E(struct params *p,const double *x,double psiR,
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

  
  tol2233b1 = fabs(no_integral_E(p,x,p->y[1][p->mpt],0)*x[1]*x[1]/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,x);
  integration_2233b1 = qromb(p->r,p->rf_fib,p->mpt,tol2233b1,&failure);
  
  if (failure) {
    *E = integration_2233b1;
    printf("failed to integrate at x = (%e,%e,%e)\n",
	   x[1],x[2],x[3]);
    return false;
  }

  *E = no_integral_E(p,x,p->y[1][p->mpt],integration_2233b1);

  return true;
}


double no_integral_E(struct params *p,const double *x,double psiR,
		     double integration_2233b1)
/*==============================================================================

  Purpose: This function computes E(x) given that the integral terms have
  already been calculated and stored in integration_2233b1.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

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
  E = 2.0/(x[1]*x[1])*integration_2233b1; 

  // add density fluctuations term
  E = (E+x[3]*x[3]*p->omega*0.5
       *(0.75*x[3]*x[3]-1));

  // add surface term tension terms
  E = E+0.5+1.0/x[1]*(-(1+p->k24)*(sin(psiR)*sin(psiR))/x[1]+2.0*p->gamma_s);  
  

  return E;
}


void compute_rf2233b1(struct params *p,const double *x)
/*==============================================================================

  Purpose: This function computes r*f_fibril(r) for all ri in r[1..mpt], and
  stores it into the array rf_fib[1..mpt].

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ============================================================================*/
{

  double f2233b1_r(struct params *p,const double *x,double ri,double sin_yi,
		   double sin_2yi,double cos_yi,double yi_p);

  int i;

  double siny, sin2y,cosy;
  
  p->rf_fib[1] = 0;
  
  for (i = 2; i <= p->mpt; i++) {  // compute f_fibril*r

    siny = sin(p->y[1][i]);

    sin2y = sin(2*p->y[1][i]);

    cosy = cos(p->y[1][i]);

    p->rf_fib[i] = p->r[i]*f2233b1_r(p,x,p->r[i],siny,sin2y,cosy,
				     p->y[2][i]);

  }

  return;
}


double f2233b1_r(struct params *p,const double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p)
/*==============================================================================

  Purpose: This function computes the integrand of the energy functional (which
  is the energy density r*f_fibril(r) in the model) at the point ri.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

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
	 *sin_yi/(ri*ri)+p->Lambda*x[3]*x[3]/4.0
	 *(4*M_PI*M_PI-x[2]*x[2]*cos_yi*cos_yi)
	 *(4*M_PI*M_PI-x[2]*x[2]*cos_yi*cos_yi));

  return ans;

}

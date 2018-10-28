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

void compute_rf2233b1(double *rf_fib,struct params *p,double *x,double *r,
		      double **y,int mpt);



// calculation of energy, given that the integral has been calculated. //

double no_integral_E(struct params *p,double *x,double psiR,
		     double integration_2233b1,int mpt);



// calculation of the energy, returns true if integral was successfully
// calculated.

bool successful_E_count(double *E,struct params *p,double *x,double *r,
			double **y,double *rf_fib,int mpt);


// files to write if calculation of E fails
void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x);
void write_QROMBfailure(double *r, double **y,double *rf_fib,int mpt,
			struct params p,double *x);
void write_SOLVDEfailure(double *r,double **y,double **y_guess,int mpt,
			 struct params p,double *x);


double f2233b1_r(struct params *p,double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp)
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

  tmp -- This is just a complicated constant passed as one variable to save
  time, tmp = 4*M_PI*M_PI/(d0*d0), where d0 is the standard d-band spacing.

  ------------------------------------------------------------------------------

  Returns: Returns r*f_fibril(r) at the point ri.

  ============================================================================*/
{
  double ans;
  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri)+p->Lambda*x[3]*x[3]/4.0
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2])
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2]));

  return ans;
}

void compute_rf2233b1(double *rf_fib,struct params *p,double *x,double *r,
		      double **y,int mpt)
/*==============================================================================

  Purpose: This function computes r*f_fibril(r) for all ri in r[1..mpt], and
  stores it into the array rf_fib[1..mpt].

  ------------------------------------------------------------------------------

  Parameters:

  rf_fib[1..mpt] -- This is the array which r*f_fibril(r) is stored in.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  r[1..mpt] -- This vector holds the grid points for psi(r).

  y[1..2][1..mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:].

  mpt -- This is the number of grid points.

  ============================================================================*/
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

double no_integral_E(struct params *p,double *x,double psiR,
		     double integration_2233b1,int mpt)
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

  mpt -- This is the number of grid points.

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



bool successful_E_count(double *E,struct params *p,double *x,double *r,
			double **y,double *rf_fib,int mpt)
/*==============================================================================

  Purpose: Given the form of psi(r) (in the array y[1..2][1..mpt]), compute the
  integrand rf_fib[1..mpt], numerically integrate this integrand, and then
  compute E(x).

  ------------------------------------------------------------------------------

  Parameters:

  *E -- pointer to variable which will store the value of E(x), if the
  calculation is successful.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  r[1..mpt] -- This vector holds the grid points for psi(r).

  y[1..2][1..mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:].

  rf_fib[1..mpt] -- This is the array which the integrand r*f_fibril is stored
  in.

  mpt -- This is the number of grid points.

  ------------------------------------------------------------------------------

  Returns: Returns true if the calculation of E(x) was successful, and stores
  E(x) in *E. Returns false if E(x) integral was not successful.

  ============================================================================*/  
  
{
  
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

  
  tol2233b1 = fabs(no_integral_E(p,x,y[1][mpt],0,mpt)*x[1]*x[1]/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(rf_fib,p,x,r,y,mpt);
  integration_2233b1 = qromb(r,rf_fib,mpt,tol2233b1,&failure);
  
  if (failure) {
    printf("failed to integrate at x = (%e,%e,%e)\n"
	   x[1],x[2],x[3],mpt);
    return false;
  }

  *E = no_integral_E(p,x,y[1][mpt],integration_2233b1,mpt);

  return true;
}

double E_calc(struct params *p,double *x,double ***c,double **s,double *r,
	      double **y,double *rf_fib,double *r_cp,double **y_cp,double conv,
	      int itmax,int *mpt,struct arr_ns *ns,int max_mpt)
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

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  c -- This is a tensor used in solvde_wrapper, c[1..2][1..2][1..max_mpt+1].

  s -- This is a tensor used by solvde_wrapper, s[1..2][1..5].

  r[1..max_mpt] -- This vector holds the grid points for psi(r). Only the 
  first mpt values, with r[mpt] = R (= x[1]) are used, but the remaining grid
  values are there in case interpolation needs to occur.

  y[1..2][1..max_mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:]. Again, only the first mpt values are used until interpolation is 
  necessary.

  rf_fib[1..max_mpt] -- This is the array which the integrand r*f_fibril is
  stored in. Again, only the first mpt values are used until interpolation is 
  necessary.

  r_cp,y_cp -- copies of r and psi(r).

  conv -- The convergence criterion for solving the ODE for psi(r) (once psi
  changes less than conv at each grid points r[1..mpt]).

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  *mpt -- Address to the integer for the number of grid points that the first
  try of calculating E(x) uses.

  *ns -- Address to structure which holds sizing of c, s, and y. This struct is
  passed to the solvde_wrapper.

  max_mpt -- The maximum number of grid points possible in r, and psi(r). If 
  interpolation is required for more than this number of grid points, then the
  calculation of E(x) is considered a failure, and the function exits to the 
  system.

  ------------------------------------------------------------------------------

  Returns: Returns the value of E(x) if the calculation is successful. If the 
  calculation is unsuccessful, the form of psi(r), dpsi/dr, and r*f_fibril(r)
  are saved in a file starting with "QROMB" and then an exit status exit(1) is
  invoked, exiting to the system.

  ============================================================================*/

{
  double h;
  double E;

  double scalv[2+1];

  int last_mpt = *mpt;

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  while ((*mpt) <= max_mpt) {

    h = x[1]/((*mpt)-1);    // compute stepsize in r[1..mpt] 

    if (need_to_interpolate(*mpt,last_mpt)) {

      printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	     x[1],x[2],x[3]);

      copy_2_arrays(r,y,r_cp,y_cp,last_mpt); // copy arrays r and y into r_cp and y_cp
      propagate_r(r,h,*mpt);
      interpolate_array(r,y,r_cp,y_cp,*mpt); // interpolate old y values (now stored in y_cp)


      last_mpt = *mpt;
    }

    else propagate_r(r,h,(*mpt));
    
    solvde_wrapper(itmax,conv,scalv,ns,*mpt,y,r,c,s,p,x,h);

    if(successful_E_count(*E,p,x,r,y,rf_fib,mpt)) return E;

    
    (*mpt) = ((*mpt)-1)*2+1;
    
  }

  // if it makes it this far, we did not successfully compute E(R)

  write_QROMBfailure(r,y,rf_fib,*mpt,*p,x); // save psi(r), rf_fib(r), and exit

  // never get here.

  return 0;
}

void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x)
/*==============================================================================

  Purpose: This function writes a string to the char array f_err. f_err will be
  the file name of a file which stores e.g. psi(r) if the calculation of E(x)
  fails.

  ------------------------------------------------------------------------------

  Parameters:

  f_err -- The char array where the file name is stored.

  err_type -- A string that is inputted to the function, and written to the
  f_err file name. The string identifies the type of error that caused the
  calculation to fail. The two possible types of errors are that the ODE
  solver solvde_wrapper fails to find psi(r) - err_type = "SOLVDE", or that
  the calculation of E(x) fails as the integration cannot be done with a 
  sufficiently low error - err_type = "QROMB".

  f_err_size -- Approximate size of the f_err char array.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Does not return a value, it just saves the file name to the f_err 
  char array.
  ============================================================================*/
  
{
  snprintf(f_err,f_err_size,"data/%s_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",err_type,
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,x[1],x[2],x[3],
	   p.gamma_s);
  return;
}

void write_QROMBfailure(double *r, double **y,double *rf_fib,int mpt,
			struct params p,double *x)
/*==============================================================================
  
  Purpose: This function saves the r, psi(r), and r*f_fibril(r) values that 
  were used in attempting to calculate E(x). This function is only called if
  the E(x) calculation fails.

  ------------------------------------------------------------------------------

  Parameters:

  r[1..max_mpt] -- This vector holds the grid points for psi(r). Only the 
  first mpt values, with r[mpt] = R (= x[1]) are used, but the remaining grid
  values are there in case interpolation needs to occur.

  y[1..2][1..max_mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:]. Again, only the first mpt values are used until interpolation is 
  necessary.

  rf_fib[1..max_mpt] -- This is the array which the integrand r*f_fibril is
  stored in. Again, only the first mpt values are used until interpolation is 
  necessary.

  *mpt -- Address to the integer for the number of grid points that the first
  try of calculating E(x) uses.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Exits with exit status exit(1);

  ============================================================================*/

{
  int i;
  FILE *broken;
  int f_err_size = 200;
  char f_err[f_err_size];

  make_f_err(f_err,"QROMB",f_err_size,p,x);

  printf("failed to integrate with qromb at x = (%e,%e,%e).\n",x[1],x[2],x[3]);
  printf("saving psi(r) shape, rf_fib(r), and exiting to system.\n");

  broken = fopen(f_err,"w");


  for (i = 1; i<=rlength; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i],rf_fib[i]);
  }
  fclose(broken);
  exit(1);
  return;
}

void write_SOLVDEfailure(double *r,double **y,double **y_guess,int mpt,
			 struct params p,double *x)
/*==============================================================================
  
  Purpose: This function saves the current forms of r and psi(r) (which has
  not relaxed to the correct psi(r) form), as well as the form used as the 
  initial guess for psi(r) which is inputted to the solvde_wrapper. This
  function is only called if solvde_wrapper fails to calculate psi(r). 

  ------------------------------------------------------------------------------

  Parameters:

  r[1..max_mpt] -- This vector holds the grid points for psi(r). Only the 
  first mpt values, with r[mpt] = R (= x[1]) are used, but the remaining grid
  values are there in case interpolation needs to occur.

  y[1..2][1..max_mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:]. Again, only the first mpt values are used until interpolation is 
  necessary.

  y_guess -- Array with the initial guess of psi(r) and dpsi/dr.

  mpt -- The number of grid points in the r and psi(r) discretization.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Exits with exit status exit(1);

  ============================================================================*/

{
  int i;
  FILE *broken1,*broken2;
  int f_err_size = 200;
  char f_err1[f_err_size];
  char f_err2[f_err_size];

  printf("failed to solve ODE when x = (%e,%e,%e).\n",x[1],x[2],x[3]);
  printf("saving current psi(r) shape (from failed solvde call), as well as the"
	 " shape of the initial guess of psi(r) (from previous call of E_calc),"
	 " and exiting to system.\n");

  make_f_err(f_err1,"SOLVDE_FAIL",f_err_size,p,x);
  broken1 = fopen(f_err1,"w");

  for (i = 1; i<=mpt; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i]);
  }
  fclose(broken1);

  make_f_err(f_err2,"SOLVDE_INITGUESS",f_err_size,p,x);
  broken2 = fopen(f_err2,"w");

  for (i = 1; i<=last_mpt; i++) {
    fprintf(broken2,"%.8e\t%.8e\t%.8e\n",r[i],y_guess[1][i],y_guess[2][i]);
  }
  fclose(broken2);

  exit(1);
  return;
}


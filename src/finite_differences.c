/* finite difference calculations for first order and second order derivatives */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]

double compute_dx(double E,double convMIN,double cushion);

double dEdxi(double E_p_dxi, double E_m_dxi,double dxi);

double ddEdxi2(double E_p_dxi, double E_m_dxi, double E, double dxi);

double ddEdxidxj(double E_p_dxi_p_dxj,double E_m_dxi_m_dxj,double E_p_dxi,
		 double E_m_dxi,double E_p_dxj,double E_m_dxj,double E,
		 double dxi, double dxj);


void derivatives_fd(double *dEdx,double E,struct params *p,double *x,
		    double ***c,double **s,double *r,double **y,double *rf_fib,
		    double *r_cp,double **y_cp,double convODE,double dx,
		    int itmax,int *mpt,struct arr_ns *ns,int max_mpt,
		    int x_size,double *hessian,bool calc_hess,double *E_p,
		    double *E_m, double *E_pij,double *E_mij)
/*==============================================================================

  Compute first order partial derivatives at the point x, and store them
  in the vector dEdx[1..x_size]. Can also compute the hessian matrix and
  store it in the flattened matrix hessian[1..x_size*x_size] if the flag
  calc_hess is set to true. Note that computing these finite-differences
  requires re-evaluating the 
  
  ============================================================================*/
{

  void compute_hessian(double *hessian,double E,struct params *p,double *x,
		       double *r,double **y,double *rf_fib,double ***c,
		       double **s,double *r_cp,double **y_cp,double convODE,
		       int itmax,int *mpt,struct arr_ns *ns,int max_mpt,
		       int x_size,double *E_p,double *E_m,double *E_pij,
		       double *E_mij,double dx);

  

  int i; 
  
  for (i = 1; i <= x_size; i++) {
    
    x[i] += dx; // add a small step
    
    E_p[i] = E_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,
		    mpt,ns,max_mpt); // compute E for this small step dx
    
    
    x[i] -= 2*dx; // remove small step, and then subtract another
    
    E_m[i] = E_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,
		    mpt,ns,max_mpt); // compute E for this small step dx
    
    x[i] += dx; // reset x to its original form
    
    dEdx[i] = (E_p[i]-E_m[i])/(2.0*dx); // compute dEdxi and store it
    
  }
  
  if (calc_hess) {
    compute_hessian(hessian,E,p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,mpt,
		    ns,max_mpt,x_size,E_p,E_m,E_pij,E_mij,dx);
  }    
  return;
  
}

void compute_hessian(double *hessian,double E,struct params *p,double *x,
		     double *r,double **y,double *rf_fib,double ***c,
		     double **s,double *r_cp,double **y_cp,double convODE,
		     int itmax,int *mpt,struct arr_ns *ns,int max_mpt,
		     int x_size,double *E_p,double *E_m,double *E_pij,
		     double *E_mij,double dx)
/*==============================================================================

  Purpose: Compute the hessian matrix of E(x) at x. Used in Newton-Raphson
  method to find minimum.

  ------------------------------------------------------------------------------

  Parameters:

  hessian[1..x_size*x_size] -- the flattened 2d hessian array, which this
  function will update.

  E -- E(x)

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

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

  c, s -- used by solvde_wrapper to find psi(r) at each xi+-dxi.

  convODE -- The convergence criterion for solving the ODE for psi(r) (once psi
  changes less than convODE at each grid points r[1..mpt]).

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  *mpt -- Address to the integer for the number of grid points that the first
  try of calculating E(x) uses.

  *ns -- Address to structure which holds sizing of c, s, and y. This struct is
  passed to the solvde_wrapper.

  max_mpt -- The maximum number of grid points possible in r, and psi(r). If 
  interpolation is required for more than this number of grid points, then the
  calculation of E(x) is considered a failure, and the function exits to the 
  system.

  x_size -- size of the array x[1..x_size].

  E_p[1..x_size], E_m[1..x_size], E_pij[1..x_size], E_mij[1..x_size] --
  dummy vectors to store E(xi+dxi), E(xi-dxi), E(xi+dxi,xj+dxj), and
  E(xi-dxi,xj-dxj) values, respectively. The only reason these matrictes are
  being passed around is because it's too computationally intensive to keep
  malloc-ing arrays each time this function is called.

  dx -- the spacing used when calculating the finite difference derivative
  approximations.

  ------------------------------------------------------------------------------

  Returns: No values returned, but the hessian matrix is calculated and stored
  in hessian[1..x_size*x_size].

  ============================================================================*/
{

  int i,j;



  for (i = 1; i <= x_size-1; i++) { // compute off-diagonal hessian terms
    for (j = i+1; j <= x_size; j++) {


      x[i] += dx; // add a small steps
      x[j] += dx;
      E_pij[i] = E_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,mpt,
			ns,max_mpt); // compute E for these small steps dx

      x[i] -= 2*dx; // remove small steps, and then subtract another small step
      x[j] -= 2*dx;

      E_mij[i] = E_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,mpt,
			ns,max_mpt); // compute E for these small steps -dx

      x[i] += dx; // reset x to its original form
      x[j] += dx;
      
      HESS(i,j) = HESS(j,i) = ddEdxidxj(E_pij[i],E_mij[i],E_p[i],E_m[i],
					E_p[j],E_m[j],E,dx,dx);
    }
  }

  for (i = 1; i <= x_size; i++) { // compute diagonal hessian terms
    HESS(i,i) = ddEdxi2(E_p[i],E_m[i],E,dx);
  }

  return;
}

double compute_dx(double E,double convMIN,double cushion)
/*==============================================================================

  Purpose: Calculates the maximum round off error by subtracting E(x+a) from
  E(x) is at most
  E(x+a)-E(x) = C*(fabs(E(x+a))+fabs(E(x))) ~ 2*C*fabs(E(x+a)) (for  small a).
  Therefore, since I'm using the usual  central difference approximation for
  derivatives, the largest error from round off in this calculation will be
  ~2*C*fabs(E)/(2*dx). That means that the size of the convergence criteria on
  gradient descent, conv, must be larger than this round off error,
  conv > 2*C*fabs(E)/(2*dx). To be safe, I've made dx = cushion*C*fabs(E)/conv,
  which ensures the error is smaller than the convergence criterion by some
  amount 1/cushion. Note that I'm not considering truncation error (from 
  finite-difference's Taylor expansion) in this calculation.

  ------------------------------------------------------------------------------

  Parameters:

  E -- E(x)

  conv -- the convergence criterion for the minimization. If all derivatives
  fabs(dEdxi) < conv, then the minimum of E(x) is said to have been reached.

  cushion -- the amount of cushion allowed between round off errors and the 
  convergence criterion (error ~ 1/cushion *convMIN).

  ------------------------------------------------------------------------------

  Returns: Returns the grid spacing dxi which will not cause significant round
  off errors.

  ============================================================================*/
  
{
  double C = 1e-16;

  return cushion*fabs(E)*C/convMIN;
}


double dEdxi(double E_p_dxi, double E_m_dxi,double dxi)
/*==============================================================================

  Purpose: Approximate partial derivative of E(x) with respect to xi by the
  centered first order finite difference with spacing dxi.

  ------------------------------------------------------------------------------

  Parameters:

  E_p_dxi -- E(xi+dxi)

  E_m_dxi -- E(xi-dxi)

  dxi -- grid spacing used to calculate finite difference.

  ------------------------------------------------------------------------------

  Returns: The approximation to dEdxi.

  ============================================================================*/
{
  return (E_p_dxi-E_m_dxi)/(2.0*dxi);
}

double ddEdxi2(double E_p_dxi, double E_m_dxi, double E, double dxi)
/*==============================================================================

  Purpose: Approximate the second order partial derivative of E(x) with respect
  to xi (along the diagonals of the Hessian matrix) by the centered second 
  order finite difference with spacing dxi.

  ------------------------------------------------------------------------------

  Parameters:

  E_p_dxi -- E(xi+dxi)

  E_m_dxi -- E(xi-dxi)

  E -- E(x)

  dxi -- grid spacing used to calculate finite difference.

  ------------------------------------------------------------------------------

  Returns: The approximation to ddEdxidxi.

  ============================================================================*/
{
  return (E_p_dxi+E_m_dxi-2*E)/(dxi*dxi);
}

double ddEdxidxj(double E_p_dxi_p_dxj,double E_m_dxi_m_dxj,double E_p_dxi,
		 double E_m_dxi,double E_p_dxj,double E_m_dxj,double E,
		 double dxi, double dxj)
/*==============================================================================

  Purpose: Approximate the second order partial derivative of E(x) (for off-
  diagonal terms in the Hessian matrix) by a second order finite difference
  that is not quite centred, but reduces the number of function evaluations of
  E(x).

  ------------------------------------------------------------------------------

  Parameters:

  E_p_dxi_p_dxj -- E(xi+dxi,xj+dxj)

  E_m_dxi_m_dxj -- E(xi-dxi,xj-dxj)

  E_p_dxi -- E(xi+dxi)

  E_m_dxi -- E(xi-dxi)

  E -- E(x)

  dxi -- grid spacing used to calculate finite difference in xi.

  dxj -- grid spacing used to calculate finite difference in xj.

  ------------------------------------------------------------------------------

  Returns: The approximation to ddEdxidxj.

  ============================================================================*/
{
  return ((E_p_dxi_p_dxj-E_p_dxi-E_p_dxj+2*E-E_m_dxi-E_m_dxj+E_m_dxi_m_dxj)
	  /(2*dxi*dxj));
}






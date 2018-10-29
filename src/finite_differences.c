/* finite difference calculations for first order and second order derivatives */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]

void derivatives_fd(double *dEdx,double E,struct params *p,double *x,
		    double ***c,double **s,double **y,double *r,double *rf_fib,
		    double **y_cp,double *r_cp,double conv,int itmax,int *mpt,
		    struct arr_ns *ns,int max_size,int x_size,double *hessian,
		    bool calc_hess,double *E_p,double *E_m, double *E_pij,
		    double *E_mij,double *x_cp)
/*==============================================================================

  Compute first order partial derivatives at the point x, and store them
  in the vector dEdx[1..x_size]. Can also compute the hessian matrix and
  store it in the flattened matrix hessian[1..x_size*x_size] if the flag
  calc_hess is set to true. Note that computing these finite-differences
  requires re-evaluating the 
  
  ============================================================================*/
{

  double E_calc(struct params *p,double *x,double *r,double **y,double *rf_fib,
		double ***c,double **s,double *r_cp,double **y_cp,double conv,
		int itmax,int *mpt,struct arr_ns *ns,int max_mpt);

  double dx;

  int i;
  
  
  
  arr_cp(x_cp,x,x_size);
  
  dx = compute_dx(E,conv);
  
  for (i = 1; i <= x_size; i++) {
    
    x_cp[i] += dx; // add a small step
    
    E_p[i] = E_calc(p,x_cp,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,
		    mpt,ns,max_size); // compute E for this small step dx
    
    
    x_cp[i] -= 2*dx; // remove small step, and then subtract another
    
    E_m[i] = E_calc(p,x_cp,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,
		    mpt,ns,max_size); // compute E for this small step dx
    
    x_cp[i] += dx; // reset x_cp to its original form
    
    dEdx[i] = (E_p[i]-E_m[i])/(2.0*dx); // compute dEdxi and store it
    
  }
  
  if (calc_hess) {
    compute_hessian(hessian,E,p,x_cp,c,s,y,r,rf_fib,y_cp,
		    r_cp,conv,itmax,mpt,last_mpt,ns,max_size,
		    x_size,E_p,E_m,E_pij,E_mij,x_cp);
  }    
  return;
  
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

double compute_dx(double E,double conv)
/*==============================================================================

  Purpose: Calculates the maximum round off error by subtracting E(x+a) from
  E(x) is at most
  E(x+a)-E(x) = C*(fabs(E(x+a))+fabs(E(x))) ~ 2*C*fabs(E(x+a)) (for  small a).
  Therefore, since I'm using the usual  central difference approximation for
  derivatives, the largest error from round off in this calculation will be
  ~2*C*fabs(E)/(2*dx). That means that the size of the convergence criteria on
  gradient descent, conv, must be larger than this round off error,
  conv > 2*C*fabs(E)/(2*dx). To be safe, I've made dx = C*fabs(E)/(0.5*conv),
  which ensures the error is half the size of the convergence criterion. Note
  that I'm not considering truncation error (from finite-difference's Taylor
  expansion) in this calculation.

  ------------------------------------------------------------------------------

  Parameters:

  E -- E(x)

  conv -- the convergence criterion for the minimization. If all derivatives
  fabs(dEdxi) < conv, then the minimum of E(x) is said to have been reached.

  ------------------------------------------------------------------------------

  Returns: Returns the grid spacing dxi which will not cause significant round
  off errors.

  ============================================================================*/
  
{
  double C = 1e-16;

  return fabs(E)*C/(0.5*conv);
}




void compute_hessian(double *hessian,double dx,double E,
		     struct params *p,double ***c,double **s,
		     double **y,double *r,double *rf_fib,
		     double **y_cp,double *r_cp,double conv,
		     int itmax,int *mpt,int last_mpt,
		     struct arr_ns *ns,int max_size,int x_size,
		     double *E_p,double *E_m, double *E_pij,
		     double *E_mij,double *x_cp)
{
  int i,j;



  for (i = 1; i <= x_size-1; i++) { // compute off-diagonal hessian terms
    for (j = i+1; j <= x_size; j++) {


      x_cp[i] += dx; // add a small steps
      x_cp[j] += dx;
      E_pij[i] = E_calc(p,x_cp,c,s,y,r,rf_fib,y_cp,r_cp,conv,itmax,
			mpt,last_mpt,ns,max_size); // compute E for these small steps dx

      x_cp[i] -= 2*dx; // remove small steps, and then subtract another small step
      x_cp[j] -= 2*dx;

      E_mij[i] = E_calc(p,x_cp,c,s,y,r,rf_fib,y_cp,r_cp,conv,itmax,
			mpt,last_mpt,ns,max_size); // compute E for these small steps -dx

      x_cp[i] += dx; // reset x_cp to its original form
      x_cp[j] += dx;
      
      HESS(i,j) = HESS(j,i) = ddEdxidxj(E_pij[i],E_mij[i],E_p[i],E_m[i],
					E_p[j],E_m[j],E,dx,dx);
    }
  }

  for (i = 1; i <= x_size; i++) { // compute diagonal hessian terms
    HESS(i,i) = ddEdxi2(E_p[i],E_m[i],E,dx);
  }

  return;
}

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"




void set_direction(double *direction,double *dEdx,double *lastdEdx,int x_size)
/*==============================================================================

  Purpose: Sets the direction of the descent and stores it in the array
  direction[1..x_size].

  ------------------------------------------------------------------------------

  Parameters:

  direction[1..x_size] -- Vector to store the direction of the descent.

  dEdx[1..x_size] -- Gradient of E(x).

  betak -- The size of the "conjugate" direction to be included in the descent
  (vs just using the gradient).

  x_size -- The size of the relevant x vector.

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, but stores the direction of the
  descent in direction[1..x_size].

  ============================================================================*/
{
  double polak_betak(double *dEdx,double *lastdEdx,int x_size);

  int i;
  double betak = polak_betak(dEdx,lastdEdx,x_size);

  for (i = 1; i <= x_size; i++) {
    direction[i] = -dEdx[i]+betak*direction[i];
  }

  return;
}

double polak_betak(double *dEdx,double *lastdEdx,int x_size)
/*==============================================================================

  Purpose: Gives the magnitude of the directional change, betak, in the
  conjugate gradient method according to Polak and Ribiere (see e.g. Numerical
  Recipes for discussion). If direction[1..x_size] is the vector for the
  direction of the next step, then direction = -dEdx+betak*last_direction.

  ------------------------------------------------------------------------------

  Parameters:

  dEdx[1..x_size] -- current value of the gradient E(x).

  lastdEdx[1..x_size] -- previous value of the gradient E(x-dx).

  ------------------------------------------------------------------------------

  Returns: betak, the magnitude of the directional change to be applied to the
  new direction.

  ============================================================================*/

{
  double vector_norm_sq(double *a,int length);

  double betak=0;
  int i;

  for (i = 1; i <= x_size; i++) betak += (dEdx[i]-lastdEdx[i])*dEdx[i];

  betak /= vector_norm_sq(lastdEdx,x_size);

  return betak;
}

double vector_norm_sq(double *a,int length)
{
  int i;
  double sum = 0;
  for (i = 1; i <= length; i++) sum += a[i]*a[i];

  return sum;
}

void armijo_backtracker(double rate,double E,double *dEdx,double *direction,
			struct params *p,double *x,double *r,double **y,
			double *rf_fib, double ***c, double **s,double *r_cp,
			double **y_cp,double conv,int itmax,int *mpt,
			struct arr_ns *ns,int max_mpt,double min_rate,
			int x_size)
/*==============================================================================

  Purpose: Given the energy, E(x), derivatives, dEdx, and descent direction,
  determine the vector x which will minimizes E.

  ------------------------------------------------------------------------------

  Parameters:

  rate -- The starting guess for the rate parameter.

  E -- E(x)

  dEdx[1..x_size] -- The gradient of E(x).

  direction[1..x_size] -- The direction of descent.

  p -- pointer to struct of all constant parameters (e.g. K33, k24).

  x[1..x_size] -- This vector holds the variable parameters x = (R,eta,delta)'.

    r[1..max_mpt] -- This vector holds the grid points for psi(r). Only the 
  first mpt values, with r[mpt] = R (= x[1]) are used, but the remaining grid
  values are there in case interpolation needs to occur.

  y[1..2][1..max_mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:]. Again, only the first mpt values are used until interpolation is 
  necessary.

  rf_fib[1..max_mpt] -- This is the array which the integrand r*f_fibril is
  stored in. Again, only the first mpt values are used until interpolation is 
  necessary.

  c -- This is a tensor used in solvde_wrapper, c[1..2][1..2][1..max_mpt+1].

  s -- This is a tensor used by solvde_wrapper, s[1..2][1..5].

  r_cp,y_cp -- copies of r and psi(r).

  conv -- The convergence criterion for solving the ODE for psi(r) (once psi
  changes less than conv at each grid points r[1..mpt]).

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  *mpt -- Address to the integer for the number of grid points that the first
  try of calculating E(x) uses.

  *ns -- Address to structure which holds sizing of c, s, and y. This struct is
  passed to the solvde_wrapper when solving the ODE for psi(r).

  max_mpt -- The maximum number of grid points possible in r, and psi(r). If 
  interpolation is required for more than this number of grid points, then the
  calculation of E(x) is considered a failure, and the function exits to the 
  system.

  min_rate -- The minimum rate allowed to be used in this calculation.

  x_size -- The number of relevant values in x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, but modifies x so the E(x) is
  minimized, given the direction of descent.

  ============================================================================*/
{

  void update_x(double *x,double rate,double *direction,int x_size);

  void reset_x(double *x,double rate,double *direction, int x_size);

  bool armijo(double E,double E_new,double rate,double *dEdx,
	      double *direction, double rho,int x_size);

  const double st = 0.7;

  const double keep_R_positive = 0.98;
  

  double E_new;
  const double rho = 0.5;

  update_x(x,rate,direction,x_size);
  
  /*
  if (x[1] <= 0 || x[1] > 10.0) {
    while ((x[1] <= 0 || x[1] > 20.0) && rate > min_rate) {
      printf("x[1] is outside of bounds, x[1] = %e\n",x[1]);
      reset_x(x,rate,direction,x_size);
      rate *= keep_R_positive;
      update_x(x,rate,direction,x_size);
    }
    return;
  }
  */
  /*
  if (x[1] <= 0) {
    reset_x(x,rate,direction,x_size);
    direction[1] = 0.0;
    update_x(x,rate,direction,x_size);
  }
  */
  E_new = F_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,mpt,ns,
		 max_mpt);

  while (!armijo(E,E_new,rate,dEdx,direction,rho,x_size)
	 && rate > min_rate) {

    reset_x(x,rate,direction,x_size);

    rate *=st;

    update_x(x,rate,direction,x_size);



    E_new = F_calc(p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,mpt,ns,
		   max_mpt);

  }
  if (rate <= min_rate) printf("rate is too small!\n");

  return;
}


void update_x(double *x,double rate,double *direction,int x_size)
/*==============================================================================
  
  Purpose: Update the x vector, given some rate and direction.

  ------------------------------------------------------------------------------

  Parameters: 

  x[1..x_size] -- The vector x = (R,eta,delta)', which is being changed to 
  minimize E(x).

  rate -- The rate at which the update to x is applied.

  direction[1..x_size] -- The magnitude and direction of change in x which
  minimizes E(x).

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, just updates x.

  ============================================================================*/
{
  int i;

  for (i = 1; i <= x_size; i++) {
    x[i] += rate*direction[i];
  }

  return;
}

void reset_x(double *x,double rate,double *direction, int x_size)
{
  int i;

  for (i =1; i <= x_size; i++) {
    x[i] -= rate*direction[i];
  }

  return;
}


  
bool armijo(double E,double E_new,double rate,double *dEdx,
	    double *direction, double rho,int x_size)
/*==============================================================================

  Purpose: Determine whether the armijo condition,
  E(x+rate*direction) <= E(x) + rho*rate*direction*dEdx, is satisfied for the
  current rate parameter.

  ------------------------------------------------------------------------------

  Parameters:

  E -- E(x)

  E_new -- E(x+rate*direction)

  rate -- The rate at which the descent updates x.

  dEdx -- The vector gradient of E(x).

  direction -- The direction of descent, as defined by the conjugate gradient
  method.

  rho -- a parameter of the armijo condition, typically 0 < rho < 1.

  x_size -- This is the number of parameters that are being varied to reach a
  minimum (i.e. the dimension of the vector dEdx).

  ------------------------------------------------------------------------------

  Returns: Returns true if the condition is satisfied, false if not.

  ============================================================================*/
{
  int i;
  double gkTdirection = 0;


  for (i = 1; i <= x_size; i++) {
    gkTdirection += dEdx[i]*direction[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdirection);


  return cond1;
}

  

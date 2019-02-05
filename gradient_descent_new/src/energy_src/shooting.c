#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define EPS 1.0e-14


bool success_brent(double *psip01,double psip02,struct params *p,
		   const double *x,double h)
/*==============================================================================

  NOTE THAT SUCCESS (TRUE) IS RETURNED EVEN IF THE ROOT IS NOT BRACKETED, as
  the form of psi(r) (stored in y) from calling this function is still reason-
  ably useful, and can be used as an initial guess for relaxation. The only time
  this function fails is if the bounds psip01 and psip02 contain the minimum,
  but it can't be found in itmax iterations.

  ============================================================================*/
{

  double F_bound(double psip0,struct params *p,const double *x,double h);

  int iter;

  int itmax = 100;

  double tol = EPS;
  double a = *psip01, b = psip02, c = b;
  double d,e,min1,min2;
  double fa = F_bound(*psip01,p,x,h);
  double fb = F_bound(psip02,p,x,h);
  double fc = fb;
  double pp,q,rr,s,tol1,psip0m;

  if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) {
    printf("root isn't bracketed!\n");
    b = *psip01;
    return true;
  }
  for (iter = 1; iter <= itmax; iter++) {
    if ((fb > 0 && fc > 0 ) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
      e = d = b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0*EPS*fabs(b)+0.5*tol;
    psip0m = 0.5*(c-b);
    if (fabs(psip0m) <= tol1 || fb == 0) {
      fb = F_bound(b,p,x,h);
      printf("psip0 = %e\n",b);
      b = *psip01;
      return true;
    }
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a==c) {
	pp = 2.0*psip0m*s;
	q = 1-s;
      } else {
	q = fa/fc;
	rr = fb/fc;
	pp = s*(2.0*psip0m*q*(q-rr)-(b-a)*(rr-1));
	q = (q-1)*(rr-1)*(s-1);
      }
      if (pp > 0) q = -q;
      pp = fabs(pp);
      min1 = 3.0*psip0m*q-fabs(tol1*q);
      min2 = fabs(e*q);
      if (2.0*pp < (min1 < min2 ? min1 : min2)) {
	e = d;
	d = pp/q;
      } else {
	d = psip0m;
	e = d;
      }
    } else {
      d = psip0m;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs(d) > tol1) b += d;
    else b += SIGN(tol1,psip0m);
    
    fb = F_bound(b,p,x,h);
  }

  printf("did not converge after %d iterations.\n",itmax);

  return false;
}



double F_bound(double psip0,struct params *p,const double *x,double h)
{

  double F_eq(double psiR,double psipR,double R,double k24);
    
  void rk4driver(struct params *p,const double *x,double h);
  
  
  double f;
  p->y[1][1] = 0;
  p->y[2][1] = psip0;
  rk4driver(p,x,h);

  f = F_eq(p->y[1][p->mpt],p->y[2][p->mpt],x[1],p->k24);

  if ((fabs(p->y[1][p->mpt])>= M_PI/2.0 || fabs(f) >=10)
      && x[3] != 0) {
    return 10;
  }
  else return f;
}


double F_solve(double psip0,struct params *p,const double *x,double h)
{
  
  double F_eq(double psiR,double psipR,double R,double k24);

  void rk4driver(struct params *p,const double *x,double h);
  
  p->y[1][1] = 0;
  p->y[2][1] = psip0;
  rk4driver(p,x,h);

  return F_eq(p->y[1][p->mpt],p->y[2][p->mpt],x[1],p->k24);
}


double F_eq(double psiR,double psipR,double R,double k24)
{
  return psipR-1-0.5*k24*sin(2*psiR)/R;
}


void rk4driver(struct params *p,const double *x,double h)
{
  
  void derivs(double *dydrj,double r_j,double *y,struct params *p,const double *x);
  void rk4(double r_j,double *y,double h,double *dydr,struct params *p,
	   const double *x,int n);

  double *ydum;
  double *dydr;
  int i;

  ydum = vector(1,2);
  dydr = vector(1,2);

  ydum[1] = p->y[1][1];
  ydum[2] = p->y[2][1];

  for (i = 2; i <= p->mpt; i++) {
    derivs(dydr,p->r[i-1],ydum,p,x);
    rk4(p->r[i-1],ydum,h,dydr,p,x,2);
    p->y[1][i] = ydum[1];
    p->y[2][i] = ydum[2];
  }

  free_vector(ydum,1,2);
  free_vector(dydr,1,2);


}



void rk4(double r_j,double *y,double h,double *dydr,struct params *p,
	 const double *x,int n)
/*==============================================================================

  Purpose: compute a fourth order runge-kutta step at the radial value r_j, and
  store the step in y. This step requires derivative information about y at r_j
  (i.e. must call derivs before using this function).

  ------------------------------------------------------------------------------

  Parameters:

  r_j -- the radial value at which dydr is computed at.

  y[1..2] -- psi and psi', i.e. y[1] = psi(r_j), y[2] = psi'(r_j).

  dydr[1..2] -- Holds the derivatives of y from the derivs call.

  p -- address pointing to the parameter struct (where K33, Lambda, etc are
  stored).

  x[1..3] -- x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, but it stores psi(r_j) and 
  psi'(r_j) in y[1] and y[2], respectively.

  ============================================================================*/
{

  void derivs(double *dydrj,double r_j,double *y,struct params *p,const double *x);
  
  double *ytmp; // stores e.g. yi+0.5*h*dyidrj
  double *dy1,*dy2;
  int i;

  ytmp = vector(1,n);
  dy1 = vector(1,n);
  dy2 = vector(1,n);

  for (i = 1; i <= n; i++) {
    ytmp[i] = y[i]+0.5*h*dydr[i];
  }

  derivs(dy1,r_j+0.5*h,ytmp,p,x);
  
  for (i = 1; i <= n; i++) {
    ytmp[i] = y[i]+0.5*h*dy1[i];
  }

  derivs(dy2,r_j+0.5*h,ytmp,p,x);

  for (i = 1; i <= n; i++) {
    ytmp[i] = y[i]+h*dy2[i];
    dy2[i] += dy1[i];         // store k2+k3, so one less array is needed.
  }

  derivs(dy1,r_j+h,ytmp,p,x);

  for (i = 1; i <= n; i++) {
    y[i] = y[i]+h/6.0*(dydr[i]+2*dy2[i]+dy1[i]);
  }

  free_vector(ytmp,1,n);
  free_vector(dy1,1,n);
  free_vector(dy2,1,n);

  return;
}

void derivs(double *dydrj,double r_j,double *y,struct params *p,const double *x)
/*==============================================================================

  Purpose: For the function dyidr = fi(r,y1,y2), where y1 is psi(r), y2 is
  psi'(r), compute the function fi(r,y1,y2) at a given r value.

  ------------------------------------------------------------------------------

  Parameters:

  dydrj[1..2] -- This vector will store the value of dydr = f(r_j,y[1],y[2]).

  r_j -- the radial value at which dydr is computed at.

  y[1..2] -- psi and psi', i.e. y[1] = psi(r_j), y[2] = psi'(r_j).

  p -- address pointing to the parameter struct (where K33, Lambda, etc are
  stored).

  x[1..3] -- x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, but stores f(r_j,y[1],y[2]) in
  the first argument dydrj.

  ============================================================================*/

{
  double siny1 = sin(y[1]);
  double sin2y1 = sin(2*y[1]);
  double cosy1 = cos(y[1]);
  double cos2y1 = cos(2*y[1]);

  if (r_j == 0) {
    dydrj[1] = y[2];
    dydrj[2] = 0;
    return;
  }

  dydrj[1] = y[2];
  dydrj[2] = (1.0/r_j*(1-y[2]-cos2y1*(1-0.5*sin2y1/r_j)
		       +p->K33*siny1*siny1*sin2y1/r_j
		       +p->Lambda*x[3]*x[3]*x[2]*x[2]*r_j
		       *(4*M_PI*M_PI-x[2]*x[2]*cosy1*cosy1)
		       *cosy1*siny1));
  return;
}

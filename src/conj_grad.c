#include <stdio.h>
#include <stdlib.h>
#include <gsl/lin_alg>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"


void update_p(struct params *p,double rate,double *arrchange,
	      double *dx, int x_size)
{
  int i;
  for (i = 1; i <= x_size; i++) dx[i] = rate*arrchange[i];

  p->R += dx[1];
  p->eta += dx[2];
  p->delta += dx[3];
  return;
}

void reset_p(struct params *p,double *dx)
{
  p->R -= dx[1];
  p->eta -= dx[2];
  p->delta -= dx[3];
  return;
}


double polak_betak(double *dEdx,double *lastdEdx)
// Sets the magnitude of the directional change from the exact gradient //
// in the conjugate gradient method according to Polak and Ribiere (see //
// e.g. Numerical Recipes for discussion). So if d_k is the new vector  //
// for the direction of the next step, d_k = -dEdx+betak*d_{k-1}. //
{
  double betak=0;
  int i;

  for (i = 1; i <= 3; i++) betak += (dEdx[i]-lastdEdx[i])*dEdx[i];

  betak /= vector_norm(lastdEdx,3);

  return betak;
}


void set_direction(double *direction,double *dEdx,double betak)
{

  direction[1] = -dEdx[1]+betak*direction[1];
  direction[2] = -dEdx[2]+betak*direction[2];
  direction[3] = -dEdx[3]+betak*direction[3];

  return;
}
  
bool armijo(double E,double E_new,double rate,double *dEdx,
	    double *direction, double rho)
{
  int i;
  double gkTdirection = 0;
  for (i = 1; i <= 3; i++) {
    gkTdirection += dEdx[i]*direction[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdirection);


  return cond1;
}

bool weak_wolfe(double E,double E_new,double rate,double *dEdx,
		double *dEdx_new,double *direction, double rho,double sigma)
{
  int i;
  double gkTdirection = 0;
  double gkTdirection_new = 0;
  for (i = 1; i <= 3; i++) {
    gkTdirection += dEdx[i]*direction[i];
    gkTdirection_new += dEdx_new[i]*direction[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdirection);
  bool cond2 = (gkTdirection_new >= sigma*gkTdirection);


  return cond1 && cond2;
}

bool strong_wolfe(double E,double E_new,double rate,double *dEdx,
		  double *dEdx_new,double *direction, double rho,double sigma)
{
  int i;
  double gkTdirection = 0;
  double gkTdirection_new = 0;
  for (i = 1; i <= 3; i++) {
    gkTdirection += dEdx[i]*direction[i];
    gkTdirection_new += dEdx_new[i]*direction[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdirection);
  bool cond2 = (fabs(gkTdirection_new) <= sigma*fabs(gkTdirection));

  //  if (cond1 && cond2) printf("strong wolfe!\n");
  return cond1 && cond2;
}

  
double armijo_backtracker(double st,double rate,double E,double *dEdx,
			  double *direction,struct params *p,double *r,double **y,
			  double *r_cp,double **y_cp,double conv,
			  int itmax,int mpt,struct arr_ns *ns,
			  int last_mpt,int max_size,double min_rate,
			  int x_size)
// Using the armijo condition (first Wolfe condition) to determine what the rate //
// parameter should be. //
{

  const double st = 0.7;
  
  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_fibc;
  double E_new;
  double *dEdx_tmp,*dx_tmp;
  double rho = 0.5;
  
  dx_tmp = vector(1,x_size);
  dEdx_tmp = vector(1,x_size);
  array_constant(0,dEdx_tmp,x_size);

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_fibc,mpt,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_mpt);

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,direction,dx_tmp,x_size);

  // calculate E, given the new values of R,eta, and delta

  // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
  single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_fibc,y_cp,r_cp,
	      dEdx_tmp,conv,itmax,&mpt,last_mpt,ns,max_size);

  last_mpt = mpt;

  reset_p(p,dx_tmp);

  while (!armijo(E,E_new,rate,dEdx,direction,rho)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,direction,dx_tmp,x_size);

    // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
    single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_fibc,y_cp,r_cp,
		dEdx_tmp,conv,itmax,&mpt,last_mpt,ns,max_size);

    last_mpt = mpt;

    reset_p(p,dx_tmp);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_fibc,mpt,ns);
  free_vector(dEdx_tmp,1,x_size);
  free_vector(dx_tmp,1,x_size);
  return rate;
}

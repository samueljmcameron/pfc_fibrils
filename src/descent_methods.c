/* Functions for computing phase field crystal model of fibrils.


 */





#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"
#include <time.h>





// conjugate gradient descent tools. //

void update_p(struct params *p,double rate,double *arrchange);

void reset_p(struct params *p,double rate,double *arrchange);

double polak_betak(double *dEdx,double *lastdEdx);

void set_dk(double *dk,double *dEdx,double betak);

bool armijo(double E,double E_new,double rate,double *dEdx,
	    double *dk, double rho);

bool weak_wolfe(double E,double E_new,double rate,double *dEdx,
		double *dEdx_new,double *dk, double rho,double sigma);

bool strong_wolfe(double E,double E_new,double rate,double *dEdx,
		  double *dEdx_new,double *dk, double rho,double sigma);

double armijo_backtracker(double st,double rate,double E,double *dEdx,
			  double *dk,struct params *p,double *r,double **y,
			  double *r_cp,double **y_cp,double conv,
			  int itmax,int npoints,struct arr_ns *ns,
			  int last_npoints,int max_size,double min_rate);

double wolfe_backtracker(double st,double rate,double E,double *dEdx,
			 double *dk,struct params *p,double *r,double **y,
			 double *r_cp,double **y_cp,double conv,
			 int itmax,int npoints,struct arr_ns *ns,
			 int last_npoints,int max_size,double min_rate,
			 bool (*wolfe)(double,double,double,double *,
				       double*,double*,double,double));

bool jumpmin(double frac_tol,double E,double *dEdR,struct params p,
	     double ***c,double **s,double **y,double *r,double *rf_,
	     double *integrand1,double *integrand2,double **y_cp,
	     double *r_cp,double conv,int itmax,int npoints,
	     int last_npoints,struct arr_ns *ns,int max_size);


void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,
	      FILE *denergyddelta,FILE *surfacetwist,
	      FILE *energydensity,double conv,int itmax,
	      int mpt,double rate0)
{
  int npoints = mpt;
  int last_npoints = mpt;
  double h;
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double *dEdx, *lastdEdx, *direction;
  double st = 0.7;
  double Elast = -1e100; // a really large, negative number to start with
  double min_rate = 1e3*conv;
  double rate = rate0;
  double lastR = 1e100;
  double lasteta = 1e100;
  double lastdelta = 1e100;
  double betak;
  int max_size = (mpt-1)*8+1;
  int count = 0;
  double frac_tol = 1e-6;
  clock_t begin = clock();
  clock_t end;

  struct arr_ns ns;
  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);
  direction = vector(1,3);

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  
  initialSlope = M_PI/(4.0*p.R);
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  printf("initial slope for guess = %lf\n", initialSlope);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  // using classical gradient descent, try to find minimum.


  array_constant(1e100,dEdx,3);
  arr_cp(lastdEdx,dEdx,3);
  array_constant(0,direction,3);



  while (non_zero_array(dEdx,conv)) {


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);

    betak = polak_betak(dEdx,lastdEdx);

    set_dk(direction,dEdx,betak);


    last_npoints = npoints;

    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.8e\t%.8e\n",p.R,dEdx[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",p.eta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",p.delta,dEdx[3]);
    fprintf(surfacetwist,"%.8e\t%.8e\n",p.R,y[1][mpt]);


    rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
			      conv,itmax,npoints,&ns,last_npoints,max_size,
			      min_rate);

    update_p(&p,rate,direction);
    count += 1;
    printf("count = %d\n",count);

    /*
    if (rate <= min_rate) {
      //      printf("last rate <= min_rate, where min_rate = %e and last "
      //     " rate = %e.\n",min_rate,rate);
      if (fabs(lastR-p.R)<conv && fabs(lasteta-p.eta)<conv
	  && fabs(lastdelta-p.delta)<conv && fabs(Elast-E)<conv) {
	//printf("changes in R, eta, and delta are smaller than %e.\n",conv);
	if(jumpmin(frac_tol,E,dEdx,p,c,s,y,r,rf_,integrand1,
		   integrand2,y_cp,r_cp,conv,itmax,npoints,
		   last_npoints,&ns,max_size)) break;
      }
    }
    */
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    Elast = E;
    


    if (count == 100000) {
      end = clock();
      printf("Elapsed after %d: %lf seconds\n",count,
	     (double)(end - begin) / CLOCKS_PER_SEC);

      exit(1);
    }

    if (p.R <= 0) {
      printf("R is being driven to negative values, R = %e!\n",p.R);
      p.R = 1e-6;
    }
      
  }
  



  printf("count = %d\n",count);
  save_psi(psi,r,y,npoints);
  save_energydensity(energydensity,r,rf_,npoints);
  printf("SAVED!\n");
  printf("E_min-E_chol = %1.2e\n",E);

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);
  free_vector(direction,1,3);
  
  return;
}


void update_p(struct params *p,double rate,double *arrchange)
{
  p->R += rate*arrchange[1];
  p->eta += rate*arrchange[2];
  p->delta += rate*arrchange[3];
  return;
}

void reset_p(struct params *p,double rate,double *arrchange)
{
  p->R -= rate*arrchange[1];
  p->eta -= rate*arrchange[2];
  p->delta -= rate*arrchange[3];
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


void set_dk(double *dk,double *dEdx,double betak)
{

  dk[1] = -dEdx[1]+betak*dk[1];
  dk[2] = -dEdx[2]+betak*dk[2];
  dk[3] = -dEdx[3]+betak*dk[3];

  return;
}
  
bool armijo(double E,double E_new,double rate,double *dEdx,
	    double *dk, double rho)
{
  int i;
  double gkTdk = 0;
  for (i = 1; i <= 3; i++) {
    gkTdk += dEdx[i]*dk[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdk);


  return cond1;
}

bool weak_wolfe(double E,double E_new,double rate,double *dEdx,
		double *dEdx_new,double *dk, double rho,double sigma)
{
  int i;
  double gkTdk = 0;
  double gkTdk_new = 0;
  for (i = 1; i <= 3; i++) {
    gkTdk += dEdx[i]*dk[i];
    gkTdk_new += dEdx_new[i]*dk[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdk);
  bool cond2 = (gkTdk_new >= sigma*gkTdk);


  return cond1 && cond2;
}

bool strong_wolfe(double E,double E_new,double rate,double *dEdx,
		  double *dEdx_new,double *dk, double rho,double sigma)
{
  int i;
  double gkTdk = 0;
  double gkTdk_new = 0;
  for (i = 1; i <= 3; i++) {
    gkTdk += dEdx[i]*dk[i];
    gkTdk_new += dEdx_new[i]*dk[i];
  }

  bool cond1 = (E_new-E <= rho*rate*gkTdk);
  bool cond2 = (fabs(gkTdk_new) <= sigma*fabs(gkTdk));

  //  if (cond1 && cond2) printf("strong wolfe!\n");
  return cond1 && cond2;
}

  
double armijo_backtracker(double st,double rate,double E,double *dEdx,
			  double *dk,struct params *p,double *r,double **y,
			  double *r_cp,double **y_cp,double conv,
			  int itmax,int npoints,struct arr_ns *ns,
			  int last_npoints,int max_size,double min_rate)
// Using the strong Wolfe conditions to determine what the rate //
// parameter should be. //
{

  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;
  double E_new;
  double *dEdx_tmp;
  double rho = 0.5;

  dEdx_tmp = vector(1,3);
  array_constant(0,dEdx_tmp,3);

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,dk);
  
  // calculate E, given the new values of R,eta, and delta

  single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,conv,itmax,&npoints,
	      last_npoints,ns,max_size);

  last_npoints = npoints;

  reset_p(p,rate,dk);

  while (!armijo(E,E_new,rate,dEdx,dk,rho)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,dk);
    single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,conv,itmax,&npoints,
		last_npoints,ns,max_size);

    last_npoints = npoints;

    reset_p(p,rate,dk);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);
  free_vector(dEdx_tmp,1,3);
  return rate;
}
  
double wolfe_backtracker(double st,double rate,double E,double *dEdx,
			 double *dk,struct params *p,double *r,double **y,
			 double *r_cp,double **y_cp,double conv,
			 int itmax,int npoints,struct arr_ns *ns,
			 int last_npoints,int max_size,double min_rate,
			 bool (*wolfe)(double,double,double,double *,
				       double*,double*,double,double))
// Using the strong Wolfe conditions to determine what the rate //
// parameter should be. //
{

  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;
  double E_new;
  double *dEdx_new;
  double rho = 0.01;
  double sigma = 0.1;

  dEdx_new = vector(1,3);
  arr_cp(dEdx_new,dEdx,3);

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,dk);
  
  // calculate E, given the new values of R,eta, and delta

  single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,conv,itmax,&npoints,
	      last_npoints,ns,max_size);

  last_npoints = npoints;

  reset_p(p,rate,dk);

  while (!wolfe(E,E_new,rate,dEdx,dEdx_new,dk,rho,sigma)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,dk);
    single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,conv,itmax,&npoints,
		last_npoints,ns,max_size);

    last_npoints = npoints;

    reset_p(p,rate,dk);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);
  free_vector(dEdx_new,1,3);
  return rate;
}

bool jumpmin(double frac_tol,double E,double *dEdx,struct params p,
	     double ***c,double **s,double **y,double *r,double *rf_,
	     double *integrand1,double *integrand2,double **y_cp,
	     double *r_cp,double conv,int itmax,int npoints,
	     int last_npoints,struct arr_ns *ns,int max_size)
{
  int i;
  double lastR = p.R;
  double lasteta = p.eta;
  double lastdelta = p.delta;
  double *tmpdEdx;
  double Elast = E;
  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;

  tmpdEdx = vector(1,3);
  arr_cp(tmpdEdx,dEdx,3);


  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);


  p.R = lastR*(1-sign(dEdx[1])*frac_tol);
  p.eta = lasteta*(1-sign(dEdx[2])*frac_tol);
  p.delta = lastdelta*(1-sign(dEdx[3])*frac_tol);

  single_calc(&E,tmpdEdx,&p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,conv,itmax,&npoints,last_npoints,
	      ns,max_size);

  printf("lastE = %e, E = %e\n",Elast,E);
  for (i = 1; i <= 3; i++) {
    printf("initial dEdx[%d] = %e, new dEdx[%d] = %e,",
	   i,dEdx[i],i,tmpdEdx[i]);
  }
  printf("\n\n");
  if (dEdx[1]*tmpdEdx[1] <=0 && dEdx[2]*tmpdEdx[2] <= 0
      && dEdx[3]*tmpdEdx[3] <= 0) {
    printf("jumped over the minimimum (derivatives all changed signs).\n");
    return true;
  }


  return false;
}


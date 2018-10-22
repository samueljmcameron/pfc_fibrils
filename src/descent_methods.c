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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>





// conjugate gradient descent tools. //

bool positive_definite(double *hessian);

void hessian_update_p(struct params *p, double *hessian, double *dEdx);

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
  double *hessian;  // note vector (vs matrix) definition of hessian
  double initialSlope;
  double E;
  double *dEdx, *lastdEdx, *direction;
  double st = 0.7;
  double lastE = -1e100; // a really large, negative number to start with
  double min_rate = conv;
  double rate = rate0;
  double lastR = p.R;
  double lasteta = p.eta;
  double lastdelta = p.delta;
  double betak;
  int max_size = (mpt-1)*8+1;
  int count = 0;
  double frac_tol = 1e-6;
  clock_t begin = clock();
  clock_t end;

  bool positivedef = false;
  int pos_def_in_a_row = 0;

  struct arr_ns ns;

  // allocate space for the arrays which will not be resized
  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);
  direction = vector(1,3);
  hessian = vector(1,3*3);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);

  // malloc the relevant arrays which may be resized
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);

  
  initialSlope = M_PI/(4.0*p.R);
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  printf("initial slope for guess = %lf\n", initialSlope);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  // using classical gradient descent, try to find minimum.


  array_constant(1e100,dEdx,3);
  arr_cp(lastdEdx,dEdx,3);
  array_constant(0,direction,3);

  
  while (pos_def_in_a_row < 20) {


    if (count % 100 != 0 || count == 0) {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		  &integrand2,y_cp,r_cp,hessian,conv,itmax,
		  &npoints,last_npoints,&ns,max_size,false);

      betak = polak_betak(dEdx,lastdEdx);

      set_dk(direction,dEdx,betak);

      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,npoints,&ns,last_npoints,max_size,
				min_rate);

      update_p(&p,rate,direction);

    } else {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		  &integrand2,y_cp,r_cp,hessian,conv,itmax,
		  &npoints,last_npoints,&ns,max_size,true);

      if (positive_definite(hessian)) {

	hessian_update_p(&p,hessian,dEdx);

	pos_def_in_a_row += 1;

	printf("Newton Raphson method worked! number of positive " 
	       "definite hessians in a row = %d!\n",
	       pos_def_in_a_row);

      } else {

	betak = polak_betak(dEdx,lastdEdx);
	
	set_dk(direction,dEdx,betak);

	rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				  conv,itmax,npoints,&ns,last_npoints,max_size,
				  min_rate);
	

	update_p(&p,rate,direction);

	pos_def_in_a_row = 0;

	printf("Newton's Raphson method didn't work. Resetting "
	       "the number of positive definite hessians in "
	       "a row back to %d.\n",pos_def_in_a_row);

      }

    }
    
    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.8e\t%.8e\n",lastR,dEdx[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.8e\t%.8e\n",lastR,y[1][mpt]);

    last_npoints = npoints;
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    lastE = E;


    count += 1;
    printf("count = %d\n",count);

    if (p.R <= 0) {
      printf("R is being driven to negative values, R = %e!\n",p.R);
      p.R = 1e-6;
    }
      
  }

  printf("\n\n\n\n\n\n"
	 "moving to (strictly) Newton-Raphson method!"
	 "\n\n\n\n\n\n");


  while (non_zero_array(dEdx,conv)) {


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		&integrand2,y_cp,r_cp,hessian,conv,itmax,
		&npoints,last_npoints,&ns,max_size,true);

    if (lastE+1e-14*fabs(lastE) < E) {      

      betak = polak_betak(dEdx,lastdEdx);
      
      set_dk(direction,dEdx,betak);

      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,npoints,&ns,last_npoints,max_size,
				min_rate);
      
      update_p(&p,rate,direction);

      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);

    } else {

      hessian_update_p(&p,hessian,dEdx);

    }


    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.8e\t%.8e\n",lastR,dEdx[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.8e\t%.8e\n",lastR,y[1][mpt]);
    

    last_npoints = npoints;
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    lastE = E;

    count += 1;
    printf("count = %d\n",count);

    if (p.R <= 0) {
      printf("R is being driven to negative values, R = %e!\n",p.R);
      p.R = 1e-6;
    }
      
  }

  end = clock();

  printf("calculation for %d iterations with Newton descent = %e\n",
	 count,(double) (end-begin)/CLOCKS_PER_SEC);


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
  free_vector(hessian,1,3*3);

  return;
}

bool positive_definite(double *hessian)
{
  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,3,3);

  gsl_matrix *mcp = gsl_matrix_alloc(3,3);
  
  gsl_matrix_memcpy(mcp,&m.matrix);

  gsl_vector *eval = gsl_vector_alloc(3);

  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(3);

  gsl_eigen_symm(mcp,eval,w);
  


  if (gsl_vector_get(eval,0) < 0 || gsl_vector_get(eval,1) < 0
      || gsl_vector_get(eval,2) < 0) {
      gsl_vector_free(eval);
      gsl_eigen_symm_free(w);
      gsl_matrix_free(mcp);
      return false;
  }
  gsl_vector_free(eval);
  gsl_eigen_symm_free(w);
  gsl_matrix_free(mcp);

  return true;
}

void hessian_update_p(struct params *p, double *hessian, double *dEdx)
{

  gsl_vector *dx = gsl_vector_alloc(3);

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,3,3);


  gsl_vector_view b = gsl_vector_view_array(dEdx+1,3);
  int s;
  gsl_permutation *perm = gsl_permutation_alloc(3);

  //  compute_eigenvalues(hessian);

  gsl_linalg_LU_decomp(&m.matrix,perm,&s);

  gsl_linalg_LU_solve(&m.matrix,perm,&b.vector,dx);

  p->R -= gsl_vector_get(dx,0);
  p->eta -= gsl_vector_get(dx,1);
  p->delta -= gsl_vector_get(dx,2);

  gsl_permutation_free(perm);
  gsl_vector_free(dx);

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

  // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
  single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,dEdx_tmp,conv,itmax,&npoints,
	      last_npoints,ns,max_size,false);

  last_npoints = npoints;

  reset_p(p,rate,dk);

  while (!armijo(E,E_new,rate,dEdx,dk,rho)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,dk);

    // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
    single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,dEdx_tmp,conv,itmax,&npoints,
		last_npoints,ns,max_size,false);

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

  // putting dEdx_new in as dummy variable for hessian, as hessian is not used
  single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,dEdx_new,conv,itmax,&npoints,
	      last_npoints,ns,max_size,false);

  last_npoints = npoints;

  reset_p(p,rate,dk);

  while (!wolfe(E,E_new,rate,dEdx,dEdx_new,dk,rho,sigma)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,dk);

    // putting dEdx_new in as dummy variable for hessian, as hessian is not used
    single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,dEdx_new,conv,itmax,&npoints,
		last_npoints,ns,max_size,false);

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
  double lastE = E;
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

  // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
  single_calc(&E,tmpdEdx,&p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,tmpdEdx,conv,itmax,&npoints,last_npoints,
	      ns,max_size,false);

  printf("lastE = %e, E = %e\n",lastE,E);
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


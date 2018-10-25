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


#define EFFECTIVE_ZERO 1e-14


// conjugate gradient descent tools. //

bool positive_eigen(gsl_vector *eigenvals, int num_x);

bool positive_definite(double *hessian, int num_x);

void hessian_update_p(struct params *p, double *hessian, double *dEdx,
		      double *dx,int num_x);

void update_p(struct params *p,double rate,double *arrchange,
	      double *dx, int num_x);


void reset_p(struct params *p,double *dx);

double polak_betak(double *dEdx,double *lastdEdx);

void set_direction(double *direction,double *dEdx,double betak);

bool armijo(double E,double E_new,double rate,double *dEdx,
	    double *direction, double rho);

bool weak_wolfe(double E,double E_new,double rate,double *dEdx,
		double *dEdx_new,double *direction, double rho,double sigma);

bool strong_wolfe(double E,double E_new,double rate,double *dEdx,
		  double *dEdx_new,double *direction, double rho,double sigma);

double armijo_backtracker(double st,double rate,double E,double *dEdx,
			  double *direction,struct params *p,double *r,double **y,
			  double *r_cp,double **y_cp,double conv,
			  int itmax,int npoints,struct arr_ns *ns,
			  int last_npoints,int max_size,double min_rate,
			  int num_x);


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
  double *dEdx, *lastdEdx, *direction, *dx;
  double st = 0.7;
  double lastE = -1e100; // a really large, negative number to start with
  double min_rate = conv;
  double rate = rate0;
  double lastR = p.R;
  double lasteta = p.eta;
  double lastdelta = p.delta;
  double betak;
  int max_size = (mpt-1)*8+1;
  int num_x0 = 3;
  int num_x = num_x0;
  int count = 0;
  double frac_tol = 1e-6;
  clock_t begin = clock();
  clock_t end;

  bool positivedef = false;
  int pos_def_in_a_row = 0;

  struct arr_ns ns;

  // allocate space for the arrays which will not be resized
  assign_ns(&ns);
  dEdx = vector(1,num_x0);
  lastdEdx = vector(1,num_x0);
  direction = vector(1,num_x0);
  dx = vector(1,num_x0);
  hessian = vector(1,num_x0*num_x0);
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


  array_constant(1e100,dEdx,num_x0);
  arr_cp(lastdEdx,dEdx,num_x0);
  array_constant(0,direction,num_x0);

  
  while (pos_def_in_a_row < 20 && fabs(dEdx[1])>EFFECTIVE_ZERO) {

    if (p.delta <= EFFECTIVE_ZERO) num_x = 1;
    else num_x = num_x0;

    if (0==1) {//(count % 100 != 0 || count == 0)) {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		  &integrand2,y_cp,r_cp,hessian,conv,itmax,
		  &npoints,last_npoints,&ns,max_size,false);

      betak = polak_betak(dEdx,lastdEdx);

      set_direction(direction,dEdx,betak);

      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,npoints,&ns,last_npoints,max_size,
				min_rate,num_x);

      update_p(&p,rate,direction,dx,num_x);

    } else {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		  &integrand2,y_cp,r_cp,hessian,conv,itmax,
		  &npoints,last_npoints,&ns,max_size,true);

      if (positive_definite(hessian,num_x)) {

	hessian_update_p(&p,hessian,dEdx,dx,num_x);

	pos_def_in_a_row += 1;

	printf("Newton Raphson method worked! number of positive " 
	       "definite hessians in a row = %d!\n",
	       pos_def_in_a_row);

      } else {

	betak = polak_betak(dEdx,lastdEdx);
	
	set_direction(direction,dEdx,betak);

	rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				  conv,itmax,npoints,&ns,last_npoints,max_size,
				  min_rate,num_x);
	

	update_p(&p,rate,direction,dx,num_x);

	pos_def_in_a_row = 0;

	printf("Newton's Raphson method didn't work. Resetting "
	       "the number of positive definite hessians in "
	       "a row back to %d.\n",pos_def_in_a_row);

      }

    }
    
    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastR,E,dEdx[1],hessian[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastR,y[1][mpt],y[2][mpt]);

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

    if (p.delta <= EFFECTIVE_ZERO) num_x = 1;
    else num_x = num_x0;

    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
		&integrand2,y_cp,r_cp,hessian,conv,itmax,
		&npoints,last_npoints,&ns,max_size,true);

    printf("hessian[1][1] = %e\n",hessian[1]);
    
    if (!positive_definite(hessian,num_x)) {      
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      betak = polak_betak(dEdx,lastdEdx);
      
      set_direction(direction,dEdx,betak);
      
      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,npoints,&ns,last_npoints,max_size,
				min_rate,num_x);
      
      update_p(&p,rate,direction,dx,num_x);
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      
    } else {
      
      hessian_update_p(&p,hessian,dEdx,dx,num_x);

    }


    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastR,E,dEdx[1],hessian[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastR,y[1][mpt],y[2][mpt]);

    last_npoints = npoints;
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    lastE = E;

    count += 1;
    printf("count = %d\n",count);


      
  }


  single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,
	      &integrand2,y_cp,r_cp,hessian,conv,itmax,
	      &npoints,last_npoints,&ns,max_size,true);

  positive_definite(hessian,num_x);

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
  free_vector(dEdx,1,num_x0);
  free_vector(lastdEdx,1,num_x0);
  free_vector(direction,1,num_x0);
  free_vector(dx,1,num_x0);
  free_vector(hessian,1,num_x0*num_x0);

  return;
}

bool positive_definite(double *hessian,int num_x)
{

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,num_x,num_x);
  gsl_matrix *mcp = gsl_matrix_alloc(num_x,num_x);
  gsl_vector *eval = gsl_vector_alloc(num_x);
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(num_x);


  gsl_matrix_memcpy(mcp,&m.matrix);
  gsl_eigen_symm(mcp,eval,w);
  
  int i;
  for (i = 0; i < num_x; i++) {
    printf("eigenvalue_%d = %.6e\n",i,gsl_vector_get(eval,i));
  }

  if (!positive_eigen(eval,num_x)) {
    
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

bool positive_eigen(gsl_vector *eigenvals, int num_x)
{
  int i;
  for (i = 0; i < num_x; i++) {
    if (gsl_vector_get(eigenvals,i) < 0) return false;
  }
  return true;
}

void hessian_update_p(struct params *p, double *hessian, double *dEdx,
		      double *dx,int num_x)
{

  gsl_vector_view dxcp = gsl_vector_view_array(dx+1,num_x);

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,num_x,num_x);

  gsl_matrix *mcp = gsl_matrix_alloc(num_x,num_x);
  
  gsl_matrix_memcpy(mcp,&m.matrix);

  gsl_vector_view b = gsl_vector_view_array(dEdx+1,num_x);
  int s;
  gsl_permutation *perm = gsl_permutation_alloc(num_x);

  //  compute_eigenvalues(hessian);

  gsl_linalg_LU_decomp(mcp,perm,&s);

  gsl_linalg_LU_solve(mcp,perm,&b.vector,&dxcp.vector);

  p->R -= gsl_vector_get(&dxcp.vector,0);
  if (p->delta > EFFECTIVE_ZERO) {
    p->eta -= gsl_vector_get(&dxcp.vector,1);
    p->delta -= gsl_vector_get(&dxcp.vector,2);
  }

  gsl_permutation_free(perm);
  gsl_matrix_free(mcp);


  return;
}


void update_p(struct params *p,double rate,double *arrchange,
	      double *dx, int num_x)
{
  int i;
  for (i = 1; i <= num_x; i++) dx[i] = rate*arrchange[i];

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
			  int itmax,int npoints,struct arr_ns *ns,
			  int last_npoints,int max_size,double min_rate,
			  int num_x)
// Using the armijo condition (first Wolfe condition) to determine what the rate //
// parameter should be. //
{

  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;
  double E_new;
  double *dEdx_tmp,*dx_tmp;
  double rho = 0.5;
  
  dx_tmp = vector(1,num_x);
  dEdx_tmp = vector(1,num_x);
  array_constant(0,dEdx_tmp,num_x);

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,direction,dx_tmp,num_x);

  // calculate E, given the new values of R,eta, and delta

  // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
  single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,dEdx_tmp,conv,itmax,&npoints,
	      last_npoints,ns,max_size,false);

  last_npoints = npoints;

  reset_p(p,dx_tmp);

  while (!armijo(E,E_new,rate,dEdx,direction,rho)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,direction,dx_tmp,num_x);

    // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
    single_calc(&E_new,dEdx_tmp,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,dEdx_tmp,conv,itmax,&npoints,
		last_npoints,ns,max_size,false);

    last_npoints = npoints;

    reset_p(p,dx_tmp);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);
  free_vector(dEdx_tmp,1,num_x);
  free_vector(dx_tmp,1,num_x);
  return rate;
}
  


/* don't need this anymore

double wolfe_backtracker(double st,double rate,double E,double *dEdx,
			 double *direction,struct params *p,double *r,double **y,
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

  dEdx_new = vector(1,num_x);
  arr_cp(dEdx_new,dEdx,num_x);

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,direction,dx,num_x);
  
  // calculate E, given the new values of R,eta, and delta

  // putting dEdx_new in as dummy variable for hessian, as hessian is not used
  single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,dEdx_new,conv,itmax,&npoints,
	      last_npoints,ns,max_size,false);

  last_npoints = npoints;

  reset_p(p,dx);

  while (!wolfe(E,E_new,rate,dEdx,dEdx_new,direction,rho,sigma)
	 && rate > min_rate) {

    rate *=st;

    update_p(p,rate,direction,dx,num_x);

    // putting dEdx_new in as dummy variable for hessian, as hessian is not used
    single_calc(&E_new,dEdx_new,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,dEdx_new,conv,itmax,&npoints,
		last_npoints,ns,max_size,false);

    last_npoints = npoints;

    reset_p(p,dx);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);
  free_vector(dEdx_new,1,num_x);
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

  tmpdEdx = vector(1,num_x);
  arr_cp(tmpdEdx,dEdx,num_x);


  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_2_arrays(r,y,rc,yc,last_npoints);


  p.R = lastR*(1-sign(dEdx[1])*frac_tol);
  p.eta = lasteta*(1-sign(dEdx[2])*frac_tol);
  p.delta = lastdelta*(1-sign(dEdx[num_x])*frac_tol);

  // putting dEdx_tmp in as dummy variable for hessian, as hessian is not used
  single_calc(&E,tmpdEdx,&p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,tmpdEdx,conv,itmax,&npoints,last_npoints,
	      ns,max_size,false);

  printf("lastE = %e, E = %e\n",lastE,E);
  for (i = 1; i <= num_x; i++) {
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

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "energy_src/nrutil.h"
#include "headerfile.h"


int delta1var_driver(double *E,struct params *p,FILE *energy)
/*==============================================================================
  This function minimizes E with respect to R, eta, and delta. It writes the
  values of R, eta, and delta which minimize E to the struct values p->R,
  p->eta, and p->delta. It also writes the energy E through the first
  argument.
  ============================================================================*/
{

  void scaledelta1var_forward(gsl_vector *y,const struct params *p);
  
  void scaledelta1var_backward(const gsl_vector *y,struct params *p);
  
  void scaledelta1var_E_backward(const double F,double *E,struct params *p);
  
  void scaledelta1var_dEdx_backward(const gsl_vector *dFdy,double *dEdx,struct params *p);

  void my_fdf(const gsl_vector *x,void *params, double *func,gsl_vector *grad);

  void df(const gsl_vector *x,void *ps,gsl_vector *g);

  double f(const gsl_vector *x_scale,void *ps);
  
  double calc_norm2(double *dEdx,int x_size);

  bool poorscaling(gsl_vector *x,int x_size);
  
  double *dEdx;

  
  dEdx = vector(1,p->x_size);

  gsl_vector *x_scale;
  x_scale = gsl_vector_alloc(p->x_size);
  scaledelta1var_forward(x_scale,p);


  
  gsl_multimin_function_fdf my_func;
  my_func.n = p->x_size;
  my_func.f = f;
  my_func.df = df;
  my_func.fdf = my_fdf;
  my_func.params = p;

  size_t iter = 0;
  size_t itermax = 100;
  int status;

  int poorscaling_count=0;
  int max_poorscaling_count = 10;

  clock_t begin = clock();



  const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc(T,p->x_size);


  gsl_multimin_fdfminimizer_set(s,&my_func,x_scale,0.001,0.01);



  printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
  	 "iteration","R","eta","delta","E","dEddelta");

  p->Escale = 0; // set initial Escale value to 0 to ensure that convergence
  //                is not obtained immediately from some large guess of Escale 

  do {
    iter ++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if (status) {
      printf ("error: %s\n", gsl_strerror (status));
      break;
    }

    status = gsl_multimin_test_gradient(s->gradient,CONV_MIN*(1+p->Escale));

    if (status == GSL_SUCCESS)
      printf("minimum found at:\n");

    scaledelta1var_backward(s->x,p);
    scaledelta1var_dEdx_backward(s->gradient,dEdx,p);
    *E = s->f;

    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,p->R,p->eta,p->delta,*E,dEdx[1]);

    if (poorscaling(s->x,p->x_size)) {
      poorscaling_count += 1;
    } else {
      poorscaling_count = 0;
    }

    if (*E>0.99*FAILED_E || poorscaling_count == max_poorscaling_count) {
      break;
    }


    p->Escale = fabs(*E);

    if (energy) {
      fprintf(energy,"%lu\t%13.6e\t%13.6e\n",iter,*E,calc_norm2(dEdx,p->x_size));
    }


  }

  while (status == GSL_CONTINUE && iter < itermax);
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x_scale);

  int returnvalue;

  if (status == GSL_SUCCESS) {
    printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
	   "iteration","R","eta","delta","E","dEddelta");

    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,p->R,p->eta,p->delta,*E,dEdx[1]);

    if (fabs(p->delta) <= DELTA_CLOSE_TO_ZERO) p->eta = sqrt(-1);

    if (p->omega == 0) p->eta=p->delta= sqrt(-1);
    
    returnvalue =  DRIVER_SUCCESS;
    
  } else {

    if (*E>0.99*FAILED_E) {

      printf("Unsuccessful calculation of energy, set to failure value E = %e\n",
	     FAILED_E);

      returnvalue = DRIVER_FAILURE;
      
    } else if (poorscaling_count == max_poorscaling_count) {

      printf("the initial guesses were not good, did not successfully find a minimum.\n");

      printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
	     "iteration","R","eta","delta","E","dEddelta");

      printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	     iter,p->R,p->eta,p->delta,*E,dEdx[1]);

      returnvalue = DRIVER_POORSCALING;

      
    } else if (iter == itermax) {

      printf("Did not successfully find a minimum. Exceeded %zu iterations.\n",iter);

      returnvalue = DRIVER_FAILURE;
      
    } else {
      
      printf("unclear why the calculation failed.\n");

      returnvalue = DRIVER_FAILURE;
      
    }
  }

  free_vector(dEdx,1,p->x_size);

  return returnvalue;
  
}

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

  for (i = 0; i < p->x_size; i++) {
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
  
  double E_calc(struct params *p);

  void scaledelta1var_backward(const gsl_vector *y,struct params *p);


  double E;
  struct params *p = ps;
  

  scaledelta1var_backward(x_scale,p);
    
  E = E_calc(p);

  return E;
}



bool poorscaling(gsl_vector *x,int x_size)
{
  int i;
  for (i = 0; i < x_size; i++) {
    if (fabs(gsl_vector_get(x,i)) > 1.0)
      return true;
  }
  return false;
}


double calc_norm2(double *dEdx,int x_size)
{
  double sum = 0;
  int i;

  for (i = 1; i <= x_size; i++) sum += dEdx[i]*dEdx[i];

  return sqrt(sum);
}

/*==============================================================================
  This re-scaling of parameters is taken from  Chapter 7.5.2 of Practical
  Optimization by Gill et al (published in 1981), on page 274, (eqn 7.17).
  ============================================================================*/


void scaledelta1var_forward(gsl_vector *y,const struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation y=D^(-1)*(x-c) from real units (x),
  to scaled values (y), where D is the transformation matrix,
  D_ii = 0.5*(xi_upper-xi_lower), and c is the constant offset vector
  c_i = 0.5*(xi_upper+xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/
{

  gsl_vector_set(y,0,2*p->delta/(p->deltaupper-p->deltalower)
		 -(p->deltaupper+p->deltalower)/(p->deltaupper-p->deltalower));
  return;
}

void scaledelta1var_backward(const gsl_vector *y,struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation x=D*y+c from scaled units (y),
  to real values (x), where D is the transformation matrix,
  D_ii = 0.5*(xi_upper-xi_lower), and c is the constant offset vector
  c_i = 0.5*(xi_upper+xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/
{

  p->delta = 0.5*((p->deltaupper-p->deltalower)*gsl_vector_get(y,0)
		  +p->deltaupper+p->deltalower);
  return;
}



void scaledelta1var_dEdx_backward(const gsl_vector *dEdy,double *dEdx,struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation dEdx = D^(-1)*dFdy from scaled
  units (dEdy) to real values (dEdx), where D is the transformation matrix,
  D = 0.5*(xi_upper-xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/

{

  dEdx[1] = 2.0/(p->deltaupper-p->deltalower)*gsl_vector_get(dEdy,0);
  return;
}

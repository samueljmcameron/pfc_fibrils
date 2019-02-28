#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "funcs/nrutil.h"
#include "headerfile.h"


int full_driver(double *E,struct params *p,FILE *energy)
/*==============================================================================
  This function minimizes E with respect to R, eta, and delta. It writes the
  values of R, eta, and delta which minimize E to the struct values p->R,
  p->eta, and p->delta. It also writes the energy E through the first
  argument.
  ============================================================================*/
{

  void scale_full_forward(gsl_vector *y,const struct params *p);
  
  void scale_full_backward(const gsl_vector *y,struct params *p);
  
  void scale_full_dEdx_backward(const gsl_vector *dFdy,double *dEdx,struct params *p);


  void my_fdf(const gsl_vector *x,void *ps, double *func,gsl_vector *grad);

  void df(const gsl_vector *x,void *ps,gsl_vector *g);

  double f(const gsl_vector *x_scale,void *ps);
  
  double calc_norm2(double *dEdx,struct params *p);

  bool poorscaling(gsl_vector *x,struct params *p);

  double *dEdx;
  
  dEdx = vector(0,p->x_size-1);

  gsl_vector *x_scale;
  x_scale = gsl_vector_alloc(p->x_size);
  scale_full_forward(x_scale,p);

  gsl_multimin_function_fdf my_func;
  my_func.n = p->x_size;
  my_func.f = f;
  my_func.df = df;
  my_func.fdf = my_fdf;
  my_func.params = p;

  size_t iter = 0;
  size_t itermax = 10000000;
  int status;

  int poorscaling_count=0;
  int max_poorscaling_count = 10;

  clock_t begin = clock();

  const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc(T,p->x_size);

  gsl_multimin_fdfminimizer_set(s,&my_func,x_scale,0.001,0.01);

  //printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
  //	 "iteration","R","eta","delta","E","dEdR","dEdeta",
  //	 "dEddelta","grad norm");

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

    scale_full_backward(s->x,p);
    scale_full_dEdx_backward(s->gradient,dEdx,p);
    *E = s->f;

    //printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
    //	   iter,p->R,p->eta,p->delta,*E,dEdx[0],dEdx[1],dEdx[2],calc_norm2(dEdx,p));

    if (poorscaling(s->x,p)) {
      poorscaling_count += 1;
    } else {
      poorscaling_count = 0;
    }

    if (*E>0.99*FAILED_E || poorscaling_count == max_poorscaling_count) {
      break;
    }


    p->Escale = fabs(*E);

    if (energy) {
      fprintf(energy,"%lu\t%13.6e\t%13.6e\n",iter,*E,calc_norm2(dEdx,p));
    }


  }

  while (status == GSL_CONTINUE && iter < itermax);
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x_scale);

  int returnvalue;

  if (status == GSL_SUCCESS) {
    printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
	   "iteration","R","eta","delta","E","dEdR","dEdeta",
	   "dEddelta","grad norm");

    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,p->R,p->eta,p->delta,*E,dEdx[0],dEdx[1],dEdx[2],calc_norm2(dEdx,p));

    if (fabs(p->delta) <= DELTA_CLOSE_TO_ZERO) p->eta = sqrt(-1);

    if (p->omega == 0) p->eta=p->delta= sqrt(-1);
    
    returnvalue =  DRIVER_SUCCESS;
    
  } else {

    if (*E>0.99*FAILED_E) {

      printf("Unsuccessful calculation of energy, set to failure value E = %e\n",
	     FAILED_E);


      printf("E = %lf\n",*E);
      printf("R = %lf\n",p->R);
      printf("eta = %lf\n",p->eta);
      printf("delta = %lf\n",p->delta);
      printf("R_c = %lf\n",p->R_c);
      printf("R_s = %lf\n",p->R_s);
      printf("psip_c = %lf\n",p->psip_c);
      printf("psip_s = %lf\n",p->psip_s);
      printf("psip_R = %lf\n",p->psip_R);
      printf("norm gradient = %lf\n",calc_norm2(dEdx,p));

      returnvalue = DRIVER_FAILURE;
      
    } else if (poorscaling_count == max_poorscaling_count) {

      printf("the initial guesses were not good, did not successfully find a minimum.\n");


      printf("E = %lf\n",*E);
      printf("R = %lf\n",p->R);
      printf("eta = %lf\n",p->eta);
      printf("delta = %lf\n",p->delta);
      printf("R_c = %lf\n",p->R_c);
      printf("R_s = %lf\n",p->R_s);
      printf("psip_c = %lf\n",p->psip_c);
      printf("psip_s = %lf\n",p->psip_s);
      printf("psip_R = %lf\n",p->psip_R);
      printf("norm gradient = %lf\n",calc_norm2(dEdx,p));

      returnvalue = DRIVER_POORSCALING;

      
    } else if (iter == itermax) {

      printf("Did not successfully find a minimum. Exceeded %zu iterations.\n",iter);

      printf("E = %lf\n",*E);
      printf("R = %lf\n",p->R);
      printf("eta = %lf\n",p->eta);
      printf("delta = %lf\n",p->delta);
      printf("R_c = %lf\n",p->R_c);
      printf("R_s = %lf\n",p->R_s);
      printf("psip_c = %lf\n",p->psip_c);
      printf("psip_s = %lf\n",p->psip_s);
      printf("psip_R = %lf\n",p->psip_R);
      printf("norm gradient = %lf\n",calc_norm2(dEdx,p));

      returnvalue = DRIVER_FAILURE;
      
    } else {
      
      printf("unclear why the calculation failed.\n");

      returnvalue = DRIVER_FAILURE;
      
    }
  }

  free_vector(dEdx,0,p->x_size-1);

  return returnvalue;
  
}


bool poorscaling(gsl_vector *x,struct params *p)
{
  int i;
  for (i = 0; i < p->x_size; i++) {
    if (fabs(gsl_vector_get(x,i)) > 1.0)
      return true;
  }
  return false;
}


double calc_norm2(double *dEdx,struct params *p)
{
  double sum = 0;
  int i;

  for (i = 0; i < p->x_size; i++) sum += dEdx[i]*dEdx[i];

  return sqrt(sum);
}


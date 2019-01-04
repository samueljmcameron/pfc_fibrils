#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "../../shared_src/nrutil.h"
#include "headerfile.h"


bool drive(double *E,struct params *p,double *x,FILE *energy)
{
  void my_fdf(const gsl_vector *x,void *params, double *func,gsl_vector *grad);

  void df(const gsl_vector *x,void *ps,gsl_vector *g);

  double f(const gsl_vector *x_scale,void *ps);
  
  double calc_norm2(double *dEdx);

  bool poorscaling(gsl_vector *x);

  double *dEdx;
  dEdx = vector(1,X_SIZE);

  gsl_vector *x_scale;
  x_scale = gsl_vector_alloc(X_SIZE);
  scale_forward(x_scale,x,p);

  gsl_multimin_function_fdf my_func;
  my_func.n = X_SIZE;
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
  s = gsl_multimin_fdfminimizer_alloc(T,X_SIZE);

  gsl_multimin_fdfminimizer_set(s,&my_func,x_scale,0.001,0.01);

  printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
	 "iteration","R","eta","delta","E","dEdR","dEdeta",
	 "dEddelta","grad norm");

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

    scale_backward(s->x,x,p);
    scale_dEdx_backward(s->gradient,dEdx,p);
    *E = s->f;
    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,x[1],x[2],x[3],*E,dEdx[1],dEdx[2],dEdx[3],calc_norm2(dEdx));
    if (poorscaling(s->x)) {
      poorscaling_count += 1;
    } else {
      poorscaling_count = 0;
    }

    if (*E>0.99*FAILED_E || poorscaling_count == max_poorscaling_count) {
      break;
    }


    p->Escale = fabs(*E);

    if (energy) {
      fprintf(energy,"%lu\t%13.6e\t%13.6e\n",iter,*E,calc_norm2(dEdx));
    }


  }

  while (status == GSL_CONTINUE && iter < itermax);
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x_scale);

  int returnvalue;

  if (status == GSL_SUCCESS) {
    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,x[1],x[2],x[3],*E,dEdx[1],dEdx[2],dEdx[3],calc_norm2(dEdx));

    if (fabs(x[3]) <= 1e-5) x[2] = sqrt(-1);

    if (p->omega == 0) x[2]=x[3]= sqrt(-1);
    
    returnvalue =  DRIVER_SUCCESS;
    
  } else {

    if (*E>0.99*FAILED_E) {

      printf("Unsuccessful calculation of energy, set to failure value E = %e\n",
	     FAILED_E);
      
    } else if (poorscaling_count == max_poorscaling_count) {

      printf("the initial guesses were not good, did not successfully find a minimum.\n");
      
    } else if (iter == itermax) {

      printf("Did not successfully find a minimum. Exceeded %zu iterations.\n",iter);
      
    } else {
      
      printf("unclear why the calculation failed.\n");
      
    }
    return false;

  }

  free_vector(dEdx,1,X_SIZE);
  
}

bool poorscaling(gsl_vector *x)
{
  int i;
  for (i = 0; i < 3; i++) {
    if (fabs(gsl_vector_get(x,i)) > 1.0)
      return true;
  }
  return false;
}


double calc_norm2(double *dEdx)
{
  double sum = 0;
  int i;

  for (i = 1; i <= X_SIZE; i++) sum += dEdx[i]*dEdx[i];

  return sqrt(sum);
}

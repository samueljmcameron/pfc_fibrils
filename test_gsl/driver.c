#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "gsl_multimin.h"
#include "nrutil.h"
#include "headerfile.h"



int main(int argc, char **argv)
{
  void my_fdf(const gsl_vector *x,void *params, double *func,gsl_vector *grad);

  void df(const gsl_vector *x,void *ps,gsl_vector *g);

  double f(const gsl_vector *x_scale,void *ps);
  
  struct params p;
  double norm;

  p.K33 = 30.0;
  p.k24 = 0.8;
  p.Lambda = 1000.0;
  p.d0 = 1.0;
  p.omega = 1000.0;
  p.Rupper = 0.04;
  p.Rlower = 0.02;
  p.etaupper = 6.3;
  p.etalower = 6.29;
  p.deltaupper = 0.83;
  p.deltalower = 0.76;
  p.gamma_s = 0.02;
  p.mpt = (MAX_M-1)/8+1;

  p.r = vector(1,MAX_M);
  p.y = matrix(1,NE,1,MAX_M);
  p.r_cp = vector(1,MAX_M);
  p.y_cp = matrix(1,NE,1,MAX_M);
  p.rf_fib = vector(1,MAX_M);
  p.s = matrix(1,NSI,1,NSJ);
  p.c = f3tensor(1,NCI,1,NCJ,1,MAX_M+1);

  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  double *x;
  gsl_vector *x_scale;
  gsl_multimin_function_fdf my_func;

  my_func.n = X_SIZE;
  my_func.f = f;
  my_func.df = df;
  my_func.fdf = my_fdf;
  my_func.params = &p;

  x_scale = gsl_vector_alloc(X_SIZE);
  x = vector(1,X_SIZE);
  x[1] = 0.03;
  x[2] = 6.295;
  x[3] = 0.815;
  scale_forward(x_scale,x,&p);

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T,X_SIZE);

  gsl_multimin_fdfminimizer_set(s,&my_func,x_scale,0.01,0.01);

  do {
    iter ++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if (status) {
      printf ("error: %s\n", gsl_strerror (status));
      break;
    }

    status = gsl_multimin_test_gradient(s->gradient,CONV_MIN);

    if (status == GSL_SUCCESS)
      printf("minimum found at:\n");

    scale_backward(s->x,x,&p);
    printf("%5lu %e %e %e %e %e %e %e %e\n",iter,x[1],x[2],x[3],s->f,
	   gsl_vector_get(s->gradient,0),
	   gsl_vector_get(s->gradient,1),
	   gsl_vector_get(s->gradient,2),
	   gsl_blas_dnrm2(s->gradient));


  }

  while (status == GSL_CONTINUE && iter < 10000);
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x_scale);

  free_vector(x,1,X_SIZE);
  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);


  return 0;
}

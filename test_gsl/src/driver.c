#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "../../shared_src/nrutil.h"
#include "headerfile.h"



int main(int argc, char **argv)
{
  void my_fdf(const gsl_vector *x,void *params, double *func,gsl_vector *grad);

  void df(const gsl_vector *x,void *ps,gsl_vector *g);

  double f(const gsl_vector *x_scale,void *ps);

  void initialize_params(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_x(double *x,struct params p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void save_psi(FILE *psi,double *r, double **y,double *rf_fib,int mpt);
  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  double *x;
  x = vector(1,X_SIZE);
  initialize_x(x,p);

  gsl_vector *x_scale;
  x_scale = gsl_vector_alloc(X_SIZE);
  scale_forward(x_scale,x,&p);

  gsl_multimin_function_fdf my_func;
  my_func.n = X_SIZE;
  my_func.f = f;
  my_func.df = df;
  my_func.fdf = my_fdf;
  my_func.params = &p;

  size_t iter = 0;
  int status;

  FILE *energy;
  initialize_file(&energy,argv[1],"energy",p);

  FILE *xvals;
  initialize_file(&xvals,argv[1],"xvals",p);

  FILE *psivsr;
  initialize_file(&psivsr,argv[1],"psivsr",p);


  clock_t begin = clock();

  const gsl_multimin_fdfminimizer_type *T;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;

  gsl_multimin_fdfminimizer *s;
  s = gsl_multimin_fdfminimizer_alloc(T,X_SIZE);

  gsl_multimin_fdfminimizer_set(s,&my_func,x_scale,0.001,0.01);

  printf("%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\t%13s\n",
	 "iteration","R","eta","delta","E","dEdR","dEdeta",
	 "dEddelta","grad norm");

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
    printf("%13lu\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	   iter,x[1],x[2],x[3],s->f,gsl_vector_get(s->gradient,0),
	   gsl_vector_get(s->gradient,1),gsl_vector_get(s->gradient,2),
	   gsl_blas_dnrm2(s->gradient));
    fprintf(energy,"%lu\t%13.6e\t%13.6e\n",iter,s->f,gsl_blas_dnrm2(s->gradient));
    fprintf(xvals,"%lu\t%13.6e\t%13.6e%13.6e\n",iter,x[1],x[2],x[3]);


  }

  while (status == GSL_CONTINUE && iter < 10000);

  clock_t end = clock();

  printf("time taken = %e\n",(double)(end-begin)/CLOCKS_PER_SEC);

  if (status == GSL_SUCCESS) {
    save_psi(psivsr,p.r,p.y,p.rf_fib,p.mpt);
  }

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

  fclose(energy);
  fclose(psivsr);
  fclose(xvals);

  return 0;
}

void initialize_x(double *x,struct params p)
{
  x[1] = p.Rguess;
  x[2] = p.etaguess;
  x[3] = p.deltaguess;
  return;
}

void initialize_param_vectors(struct params *p)
{

  void linearGuess(double *r, double **y, double initialSlope,double h,int mpt);

  p->r = vector(1,MAX_M);
  p->y = matrix(1,NE,1,MAX_M);
  p->r_cp = vector(1,MAX_M);
  p->y_cp = matrix(1,NE,1,MAX_M);
  p->rf_fib = vector(1,MAX_M);
  p->s = matrix(1,NSI,1,NSJ);
  p->c = f3tensor(1,NCI,1,NCJ,1,MAX_M+1);


  double h = p->Rguess/(p->mpt-1);
  double slopeguess = M_PI/(4.0*p->Rguess);
  linearGuess(p->r,p->y,slopeguess,h,p->mpt); //linear initial guess for psi(r)

  return;
}

void initialize_params(struct params *p,char **args)
{


  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->d0);
  sscanf(args[6],"%lf",&p->omega);
  sscanf(args[7],"%lf",&p->gamma_s);
  sscanf(args[8],"%lf",&p->Rguess);
  sscanf(args[9],"%lf",&p->etaguess);
  sscanf(args[10],"%lf",&p->deltaguess);
  sscanf(args[11],"%lf",&p->Rupper);
  sscanf(args[12],"%lf",&p->Rlower);
  sscanf(args[13],"%lf",&p->etaupper);
  sscanf(args[14],"%lf",&p->etalower);
  sscanf(args[15],"%lf",&p->deltaupper);
  sscanf(args[16],"%lf",&p->deltalower);

  p->mpt = (MAX_M-1)/8+1;

  printf("parameter values for minimization:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
  printf("d0 = %e\n",p->d0);
  printf("omega = %e\n",p->omega);
  printf("gamma_s = %e\n",p->gamma_s);


  printf("initial guesses for R, eta, and delta:\n");

  printf("guess for R = %e\n",p->Rguess);
  printf("guess for eta = %e\n",p->etaguess);
  printf("guess for delta = %e\n",p->deltaguess);


  printf("estimated ranges of R, eta, and delta:\n");

  printf("upper bound estimate of R = %e\n",p->Rupper);
  printf("lower bound estimate of R = %e\n",p->Rlower);
  printf("upper bound estimate of eta = %e\n",p->etaupper);
  printf("lower bound estimate of eta = %e\n",p->etalower);
  printf("upper bound estimate of delta = %e\n",p->deltaupper);
  printf("lower bound estimate of delta = %e\n",p->deltalower);

  return;

}

void initialize_file(FILE **output,char *path,char *fname,struct params p)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e.txt",p.K33,p.k24,p.Lambda,p.d0,p.omega,p.gamma_s);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}

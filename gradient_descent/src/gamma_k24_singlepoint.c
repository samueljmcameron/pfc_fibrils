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


int main(int argc, char **argv)
{
  void initialize_params(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_x(double *x,struct params p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void save_psivsr(FILE *psivsr,struct params *p);

  void set_x_NAN(double *E,double *x,int xsize);

  void save_observables(FILE *observables,double E,double *x,struct params *p);

  void reset_guess_vals(double *x, struct params *p);
  
  int drive(double *E,struct params *p,double *x,FILE *energy);
  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  double *x;
  x = vector(1,X_SIZE);
  initialize_x(x,p);

  FILE *observables;
  initialize_file(&observables,argv[1],"observables",p);

  FILE *psivsr;
  initialize_file(&psivsr,argv[1],"psivsr",p);

  double E;

  int calculation = drive(&E,&p,x,(NULL));

  
  if (calculation == DRIVER_POORSCALING) {
    printf("RETRYING!\n");
    reset_guess_vals(x,&p);
    calculation = drive(&E,&p,x,(NULL));
  }
  if (calculation == DRIVER_SUCCESS) {
    
    save_psivsr(psivsr,&p);
    
  } else {
    
    set_x_NAN(&E,x,X_SIZE);
    
  }

  save_observables(observables,E,x,&p);

  free_vector(x,1,X_SIZE);
  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(psivsr);
  fclose(observables);

  return 0;

}


void reset_guess_vals(double *x, struct params *p)
/*==============================================================================

  Purpose: This function is called if the scaling of x is poor (i.e. if scaledx
  values in the minimization procedure are consistently not between -1 and 1 for
  at least one of the x components (usually it is x[1]). The function assumes
  that the current values of x are okay (as they've been optimized somewhat by
  the optimization procedure) and only the scaling itself needs to be reset.
  ============================================================================*/
{
  p->Rguess = x[1];
  p->Rupper = 1.5*p->Rguess;
  p->Rlower = 0.75*p->Rguess;

  if (x[2] == x[2]) {
    p->etaguess = x[2];
    p->etaupper = p->etaguess+0.1;
    p->etalower = p->etaguess-0.02;
  }
  if (x[3] == x[3] && fabs(x[3])>DELTA_CLOSE_TO_ZERO) {
    p->deltaguess = x[3];
    p->deltaupper = 0.818;
    if (p->deltaguess < 0.81) {
      p->deltalower = 0.95*p->deltaguess;
    } else {
      p->deltalower = 0.81;
    }
  }
  
  return;
}

void save_psivsr(FILE *psivsr,struct params *p)
{
  int i;

  for (i = 1; i <= p->mpt; i++) {
    fprintf(psivsr,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",p->r[i],p->y[1][i],p->y[2][i],
	    p->rf_fib[i]);
  }
  printf("psi(R) = %1.2e\n",p->y[1][p->mpt]);
  return;
}

void save_observables(FILE *observables,double E,double *x,struct params *p)
{

  fprintf(observables,"%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	  E,x[1],x[2],x[3],p->y[1][p->mpt]);
  return;
}


void set_x_NAN(double *E,double *x,int xsize)
{
  int i;
  for (i = 1; i <= xsize; i++) x[i] = sqrt(-1);
  *E = sqrt(-1);
  return;
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

  printf("\n\nparameter values for minimization:\n");

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

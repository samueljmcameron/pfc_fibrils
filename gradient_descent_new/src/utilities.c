#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "headerfile.h"
#include "energy_src/nrutil.h"


void reset_guess_vals(struct params *p)
/*==============================================================================

  Purpose: This function is called if the scaling of R,eta,delta is poor (i.e.
  if scaledx values in the minimization procedure are consistently not between 
  -1 and 1 for at least one of the R,eta,delta components (usually it is R).
  The function assumes that the current values of R,eta,delta are okay (as 
  they've been optimized somewhat by the optimization procedure) and only the
  scaling itself needs to be reset.
  ============================================================================*/
{

  double current_R = p->R;
  double current_eta = p->eta;
  double current_delta = p->delta;
  
  p->Rguess = current_R;
  p->Rupper = 1.5*p->Rguess;
  p->Rlower = 0.75*p->Rguess;

  if (current_eta == current_eta) { // if eta is not NAN
    p->etaguess = current_eta;
    p->etaupper = p->etaguess+0.1;
    p->etalower = p->etaguess-0.02;
  }
  if (current_delta == current_delta && fabs(current_delta)>DELTA_CLOSE_TO_ZERO) {
    // if delta is not NAN and not close to zero
    p->deltaguess = current_delta;
    p->deltaupper = 0.818;
    if (p->deltaguess < 0.81) {
      p->deltalower = 0.95*p->deltaguess;
    } else {
      p->deltalower = 0.81;
    }
  }
  
  return;
}


void save_observables(FILE *observables,double E,struct params *p)
{

  fprintf(observables,"%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	  E,p->R,p->eta,p->delta,p->y[1][p->mpt]);
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



void set_NAN(double *E,struct params *p)
{

  p->R = sqrt(-1);
  p->eta = sqrt(-1);
  p->delta = sqrt(-1);
  *E = sqrt(-1);
  
  return;
}

void initialize_R_eta_delta(struct params *p)
{
  p->R = p->Rguess;
  p->eta = p->etaguess;
  p->delta = p->deltaguess;
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
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->Rguess);
  sscanf(args[8],"%lf",&p->etaguess);
  sscanf(args[9],"%lf",&p->deltaguess);
  sscanf(args[10],"%lf",&p->Rupper);
  sscanf(args[11],"%lf",&p->Rlower);
  sscanf(args[12],"%lf",&p->etaupper);
  sscanf(args[13],"%lf",&p->etalower);
  sscanf(args[14],"%lf",&p->deltaupper);
  sscanf(args[15],"%lf",&p->deltalower);

  p->mpt = (MAX_M-1)/8+1;

  printf("\n\nparameter values for minimization:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
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

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e.txt",p.K33,p.k24,p.Lambda,p.omega,p.gamma_s);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}

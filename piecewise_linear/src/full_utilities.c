#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "headerfile.h"
#include "funcs/nrutil.h"


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
  
  p->Rguess = p->R;
  p->Rupper = 1.5*p->Rguess;
  p->Rlower = 0.75*p->Rguess;

  p->R_cguess = p->R_c;
  p->R_cupper = 1.5*p->R_cguess;
  p->R_clower = 0.75*p->R_cguess;

  p->R_sguess = p->R_s;
  p->R_supper = 1.5*p->R_sguess;
  p->R_slower = 0.75*p->R_sguess;

  p->psip_cguess = p->psip_c;
  p->psip_cupper = 1.5*p->psip_cguess;
  p->psip_clower = 0.75*p->psip_cguess;

  p->psip_sguess = p->psip_s;
  p->psip_supper = 1.5*p->psip_sguess;
  p->psip_slower = 0.75*p->psip_sguess;

  p->psip_Rguess = p->psip_R;
  p->psip_Rupper = 1.5*p->psip_Rguess;
  p->psip_Rlower = 0.75*p->psip_Rguess;


  if (p->eta == p->eta) { // if eta is not NAN
    p->etaguess = p->eta;
    p->etaupper = p->etaguess+0.1;
    p->etalower = p->etaguess-0.02;
  }
  if (p->delta == p->delta && fabs(p->delta)>DELTA_CLOSE_TO_ZERO) {
    // if delta is not NAN and not close to zero
    p->deltaguess = p->delta;
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

  fprintf(observables,"%13.6e\t%13.6e\t%13.6e\t%13.6e\t%13.6e\t"
	  "%13.6e\t%13.6e\t%13.6e\t%13.6e\n",
	  E,p->R,p->eta,p->delta,p->R_c,p->R_s,p->psip_c,
	  p->psip_s,p->psip_R);
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

void initialize_x(struct params *p)
{
  p->R = p->Rguess;
  p->eta = p->etaguess;
  p->delta = p->deltaguess;
  p->R_c = p->R_cguess;
  p->R_s = p->R_sguess;
  p->psip_c = p->psip_cguess;
  p->psip_s = p->psip_sguess;
  p->psip_R = p->psip_Rguess;
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
  sscanf(args[10],"%lf",&p->R_cguess);
  sscanf(args[11],"%lf",&p->R_sguess);
  sscanf(args[12],"%lf",&p->psip_cguess);
  sscanf(args[13],"%lf",&p->psip_sguess);
  sscanf(args[14],"%lf",&p->psip_Rguess);
  sscanf(args[15],"%lf",&p->Rupper);
  sscanf(args[16],"%lf",&p->Rlower);
  sscanf(args[17],"%lf",&p->etaupper);
  sscanf(args[18],"%lf",&p->etalower);
  sscanf(args[19],"%lf",&p->deltaupper);
  sscanf(args[20],"%lf",&p->deltalower);
  sscanf(args[21],"%lf",&p->R_cupper);
  sscanf(args[22],"%lf",&p->R_clower);
  sscanf(args[23],"%lf",&p->R_supper);
  sscanf(args[24],"%lf",&p->R_slower);
  sscanf(args[25],"%lf",&p->psip_cupper);
  sscanf(args[26],"%lf",&p->psip_clower);
  sscanf(args[27],"%lf",&p->psip_supper);
  sscanf(args[28],"%lf",&p->psip_slower);
  sscanf(args[29],"%lf",&p->psip_Rupper);
  sscanf(args[30],"%lf",&p->psip_Rlower);

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
  printf("upper bound estimate of R_c = %e\n",p->R_cupper);
  printf("lower bound estimate of R_c = %e\n",p->R_clower);
  printf("upper bound estimate of R_s = %e\n",p->R_supper);
  printf("lower bound estimate of R_s = %e\n",p->R_slower);
  printf("upper bound estimate of psip_c = %e\n",p->psip_cupper);
  printf("lower bound estimate of psip_c = %e\n",p->psip_clower);
  printf("upper bound estimate of psip_s = %e\n",p->psip_supper);
  printf("lower bound estimate of psip_s = %e\n",p->psip_slower);
  printf("upper bound estimate of psip_R = %e\n",p->psip_Rupper);
  printf("lower bound estimate of psip_R = %e\n",p->psip_Rlower);

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

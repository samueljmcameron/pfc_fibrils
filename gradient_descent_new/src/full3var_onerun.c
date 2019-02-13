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

  void initialize_R_eta_delta(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void save_psivsr(FILE *psivsr,struct params *p);

  void set_NAN(double *E,struct params *p);

  void save_observables(FILE *observables,double E,struct params *p);

  void reset_guess_vals(struct params *p);
  
  int drive(double *E,struct params *p,FILE *energy);
  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);


  initialize_R_eta_delta(&p);

  FILE *observables;
  initialize_file(&observables,argv[1],"observables",p);


  double E;

  int calculation = drive(&E,&p,(NULL));

  
  if (calculation == DRIVER_POORSCALING) {
    printf("RETRYING!\n");
    reset_guess_vals(&p);
    calculation = drive(&E,&p,(NULL));
  }
  if (calculation == DRIVER_SUCCESS) {
    printf("success!\n");
    
  } else {
    
    set_NAN(&E,&p);
    
  }

  save_observables(observables,E,&p);


  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(observables);

  return 0;

}


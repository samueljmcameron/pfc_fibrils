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


int main(int argc, char **argv)
{
  void initialize_params(struct params *p,char **args);

  void initialize_x(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void save_psivsr(FILE *psivsr,struct params *p);

  void set_NAN(double *E,struct params *p);

  void save_observables(FILE *observables,double E,struct params *p);

  void reset_guess_vals(struct params *p);
  
  int full_driver(double *E,struct params *p,FILE *energy);

  printf("hello\n");
  
  struct params p; 
  initialize_params(&p,argv);

  printf("hello\n");
  
  p.x_size = 8;

  initialize_x(&p);

  FILE *observables;
  initialize_file(&observables,argv[1],"observables",p);


  double E;

  int calculation = full_driver(&E,&p,(NULL));

  
  if (calculation == DRIVER_POORSCALING) {
    printf("RETRYING!\n");
    reset_guess_vals(&p);
    calculation = full_driver(&E,&p,(NULL));
  }
  if (calculation == DRIVER_SUCCESS) {
    printf("success!\n");
    
  } else {
    
    set_NAN(&E,&p);
    
  }

  save_observables(observables,E,&p);

  fclose(observables);

  return 0;

}


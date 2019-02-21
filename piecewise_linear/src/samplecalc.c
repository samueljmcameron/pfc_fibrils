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

  void setparams(struct params *p);

  double Efunc(struct params *p);

  struct params p;

  setparams(&p);

  printf("E = %lf\n",Efunc(&p));
  
  return 0;

}

void setparams(struct params *p)
{

  p->K33 = 30.0;
  p->gamma_s = 0.04;
  p->k24 = 0.1;
  p->Lambda = 1000.0;
  p->omega = 5.0;
  
  p->R = 0.88;
  p->eta = 6.32;
  p->delta = 0.813;

  p->R_c = 0.05;
  p->R_s = 0.865;
  p->psip_c = 2.0;
  p->psip_s = 0.0;
  p->psip_R = 0.33;

  
  return;
}

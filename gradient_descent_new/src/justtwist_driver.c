#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "energy_src/nrutil.h"
#include "headerfile.h"


double justtwist_driver(struct params *p,FILE *energy)
/*==============================================================================
  This function minimizes E with respect to psi(r) for fixed R, eta, and
  delta values.
  ============================================================================*/
{

  double E_calc(struct params *p);

  return E_calc(p);
  
}

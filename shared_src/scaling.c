#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "headerfile.h"


void scale_forward(gsl_vector *y,const double *x,struct params *p)
{
  gsl_vector_set(y,0,2*x[1]/(p->Rupper-p->Rlower)
		 -(p->Rupper+p->Rlower)/(p->Rupper-p->Rlower));
  gsl_vector_set(y,1,2*x[2]/(p->etaupper-p->etalower)
		 -(p->etaupper+p->etalower)/(p->etaupper-p->etalower));
  gsl_vector_set(y,2,2*x[3]/(p->deltaupper-p->deltalower)
		 -(p->deltaupper+p->deltalower)/(p->deltaupper-p->deltalower));
  return;
}

void scale_backward(const gsl_vector *y, double *x,struct params *p)
{
  x[1] = 0.5*((p->Rupper-p->Rlower)*gsl_vector_get(y,0)
	      +p->Rupper+p->Rlower);
  x[2] = 0.5*((p->etaupper-p->etalower)*gsl_vector_get(y,1)
	      +p->etaupper+p->etalower);
  x[3] = 0.5*((p->deltaupper-p->deltalower)*gsl_vector_get(y,2)
	      +p->deltaupper+p->deltalower);
  return;
}

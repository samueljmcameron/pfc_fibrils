/*==============================================================================
  This re-scaling of parameters is taken from  Chapter 7.5.2 of Practical
  Optimization by Gill et al (published in 1981), on page 274, (eqn 7.17).
  ============================================================================*/


#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "headerfile.h"


void scale_forward(gsl_vector *y,const struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation y=D^(-1)*(x-c) from real units (x),
  to scaled values (y), where D is the transformation matrix,
  D_ii = 0.5*(xi_upper-xi_lower), and c is the constant offset vector
  c_i = 0.5*(xi_upper+xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/
{
  gsl_vector_set(y,0,2*p->R/(p->Rupper-p->Rlower)
		 -(p->Rupper+p->Rlower)/(p->Rupper-p->Rlower));
  gsl_vector_set(y,1,2*p->eta/(p->etaupper-p->etalower)
		 -(p->etaupper+p->etalower)/(p->etaupper-p->etalower));
  gsl_vector_set(y,2,2*p->delta/(p->deltaupper-p->deltalower)
		 -(p->deltaupper+p->deltalower)/(p->deltaupper-p->deltalower));
  return;
}

void scale_backward(const gsl_vector *y,struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation x=D*y+c from scaled units (y),
  to real values (x), where D is the transformation matrix,
  D_ii = 0.5*(xi_upper-xi_lower), and c is the constant offset vector
  c_i = 0.5*(xi_upper+xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/
{
  p->R = 0.5*((p->Rupper-p->Rlower)*gsl_vector_get(y,0)
	      +p->Rupper+p->Rlower);
  p->eta = 0.5*((p->etaupper-p->etalower)*gsl_vector_get(y,1)
		+p->etaupper+p->etalower);
  p->delta = 0.5*((p->deltaupper-p->deltalower)*gsl_vector_get(y,2)
		  +p->deltaupper+p->deltalower);
  return;
}



void scale_dEdx_backward(const gsl_vector *dEdy,double *dEdx,struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation dEdx = D^(-1)*dFdy from scaled
  units (dEdy) to real values (dEdx), where D is the transformation matrix,
  D = 0.5*(xi_upper-xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/

{
  dEdx[1] = 2.0/(p->Rupper-p->Rlower)*gsl_vector_get(dEdy,0);
  dEdx[2] = 2.0/(p->etaupper-p->etalower)*gsl_vector_get(dEdy,1);
  dEdx[3] = 2.0/(p->deltaupper-p->deltalower)*gsl_vector_get(dEdy,2);
  return;
}

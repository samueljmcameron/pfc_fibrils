#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "headerfile.h"

void scale_full_forward(gsl_vector *y,const struct params *p)
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
  gsl_vector_set(y,3,2*p->R_c/(p->R_cupper-p->R_clower)
		 -(p->R_cupper+p->R_clower)/(p->R_cupper-p->R_clower));
  gsl_vector_set(y,4,2*p->R_s/(p->R_supper-p->R_slower)
		 -(p->R_supper+p->R_slower)/(p->R_supper-p->R_slower));
  gsl_vector_set(y,5,2*p->psip_c/(p->psip_cupper-p->psip_clower)
		 -(p->psip_cupper+p->psip_clower)/(p->psip_cupper-p->psip_clower));
  gsl_vector_set(y,6,2*p->psip_s/(p->psip_supper-p->psip_slower)
		 -(p->psip_supper+p->psip_slower)/(p->psip_supper-p->psip_slower));
  gsl_vector_set(y,7,2*p->psip_R/(p->psip_Rupper-p->psip_Rlower)
		 -(p->psip_Rupper+p->psip_Rlower)/(p->psip_Rupper-p->psip_Rlower));
  return;
}

void scale_full_backward(const gsl_vector *y, struct params *p)
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
  p->R_c = 0.5*((p->R_cupper-p->R_clower)*gsl_vector_get(y,3)
		  +p->R_cupper+p->R_clower);
  p->R_s = 0.5*((p->R_supper-p->R_slower)*gsl_vector_get(y,4)
		  +p->R_supper+p->R_slower);
  p->psip_c = 0.5*((p->psip_cupper-p->psip_clower)*gsl_vector_get(y,5)
		  +p->psip_cupper+p->psip_clower);
  p->psip_s = 0.5*((p->psip_supper-p->psip_slower)*gsl_vector_get(y,6)
		  +p->psip_supper+p->psip_slower);
  p->psip_R = 0.5*((p->psip_Rupper-p->psip_Rlower)*gsl_vector_get(y,7)
		  +p->psip_Rupper+p->psip_Rlower);
  return;
}

void scale_full_dEdx_backward(const gsl_vector *dEdy,double *dEdx,struct params *p)
/*==============================================================================
  Purpose: Make the scaling transformation dEdx = D^(-1)*dFdy from scaled
  units (dEdy) to real values (dEdx), where D is the transformation matrix,
  D = 0.5*(xi_upper-xi_lower), where xi_upper is e.g. p->Rupper if i=1.

  ============================================================================*/

{
  dEdx[0] = 2.0/(p->Rupper-p->Rlower)*gsl_vector_get(dEdy,0);
  dEdx[1] = 2.0/(p->etaupper-p->etalower)*gsl_vector_get(dEdy,1);
  dEdx[2] = 2.0/(p->deltaupper-p->deltalower)*gsl_vector_get(dEdy,2);
  dEdx[3] = 2.0/(p->R_cupper-p->R_clower)*gsl_vector_get(dEdy,3);
  dEdx[4] = 2.0/(p->R_supper-p->R_slower)*gsl_vector_get(dEdy,4);
  dEdx[5] = 2.0/(p->psip_cupper-p->psip_clower)*gsl_vector_get(dEdy,5);
  dEdx[6] = 2.0/(p->psip_supper-p->psip_slower)*gsl_vector_get(dEdy,6);
  dEdx[7] = 2.0/(p->psip_Rupper-p->psip_Rlower)*gsl_vector_get(dEdy,7);
  
  return;
}


void my_fdf(const gsl_vector *x_scale, void *ps, double *func, gsl_vector *grad)
{
  double f(const gsl_vector *x_scale,void *ps);
  void df(const gsl_vector *x_scale,void *ps,gsl_vector *g);

  *func = f(x_scale,ps);
  df(x_scale,ps,grad);
  return;
  
}

void df(const gsl_vector *x_scale,void *ps,gsl_vector *g)
{

  void scale_full_backward(const gsl_vector *y, struct params *p);

  double dEdR(struct params *p);
  
  double dEdeta(struct params *p);

  double dEddelta(struct params *p);

  double dEdR_c(struct params *p);
  
  double dEdR_s(struct params *p);

  double dEdpsip_c(struct params *p);

  double dEdpsip_s(struct params *p);

  double dEdpsip_R(struct params *p);

  struct params *p = ps;

  scale_full_backward(x_scale,p);

  gsl_vector_set(g,0,dEdR(p));

  gsl_vector_set(g,1,dEdeta(p));

  gsl_vector_set(g,2,dEddelta(p));

  gsl_vector_set(g,3,dEdR_c(p));

  gsl_vector_set(g,4,dEdR_s(p));

  gsl_vector_set(g,5,dEdpsip_c(p));

  gsl_vector_set(g,6,dEdpsip_s(p));

  gsl_vector_set(g,7,dEdpsip_R(p));

  return;

}

double f(const gsl_vector *x_scale,void *ps)
{

  double Efunc(struct params *p);
  
  void scale_full_backward(const gsl_vector *y, struct params *p);

  
  double E;

  struct params *p = ps;
  


  scale_full_backward(x_scale,p);

  
  E = Efunc(p);

  return E;
  
}

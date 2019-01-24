#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "headerfile.h"



double ufunc(double x_1,double x_2,double zeta)
{
  return (1-zeta)*(1-zeta)*(x_2*x_2-x_1*x_1);
}

double vfunc(double x_1,double x_2,double xi,double zeta)
{
 
  double ans;
  
  if (fabs(zeta)< SMALL) {

    double a0,a1,a2,a3;

    a0 = -2*sin(2*xi)*(x_2-x_1);

    a1 = 2*sin(2*xi)*(x_2-x_1)-2*cos(2*xi)*(x_2*x_2-x_1*x_1);

    a2 = 2*cos(2*xi)*(x_2*x_2-x_1*x_1)+4.0/3.0*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1);

    a3 = (-4.0/3.0*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1)
	  +2.0/3.0*np.cos(2*xi)*(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1));
    
    ans = a0+a1*zeta+a2*zeta*zeta+a3*zeta*zeta*zeta;
    
  } else {

    ans = cos(2*(zeta*x_2+xi))-cos(2*(zeta*x_1+xi));

    ans *= (1-zeta)/zeta;

  }

  return ans;

}

double f_1func(double x_1,double x_2, double xi,double zeta)
{

  double f_1_integrand(double u,void *params);

  double ans;

  if (fabs(zeta)<SMALL) {

    double a0,a1,a2,a3;

    if ((x_1 < ZERO || x_2 < ZERO) && xi < ZERO) {
      
      a0 = 0.0;

    } else {
      
      a0 = sin(2*xi)*sin(2*xi)*log(x_2/x_1);

    }

    a1 = 4*(x_2-x_1)*cos(2*xi)*sin(2*xi);

    a2 = 2*(x_2*x_2-x_1*x_1)*(cos(2*xi)*cos(2*xi)-sin(2*xi)*sin(2*xi));

    a3 = -32.0/9.0*(x_2*x_2*x_2-x_1*x_1*x_1)*sin(2*xi)*cos(2*xi);

    ans = a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta;

  } else {

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

    double error;

    gsl_function F;

    double *v = malloc(2*sizeof(double));

    v[0] = xi;
    v[1] = zeta;

    F.function = &f_1_integrand;
    F.params = v;

    gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);

    gsl_integration_workspace_free (w);

    free(v);

  }

  return ans;

}

double f_2func(double x_1,double x_2, double xi,double zeta)
{

  double f_2_integrand(double u,void *params);

  double ans;

  if (fabs(zeta)<SMALL) {

    double a0,a1,a2,a3;

    if ((x_1 < ZERO || x_2 < ZERO) && xi < ZERO) {
      
      a0 = 0.0;

    } else {
      
      a0 = sin(xi)*sin(xi)*sin(xi)*sin(xi)*log(x_2/x_1);

    }

    a1 = 4*(x_2-x_1)*cos(xi)*sin(xi)*sin(xi)*sin(xi);

    a2 = (x_2*x_2-x_1*x_1)*sin(xi)*sin(xi)*(3*cos(xi)*cos(xi)-sin(xi)*sin(xi));

    a3 = (4.0/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*sin(xi)*cos(xi)
	  *(cos(xi)*cos(xi)-5*sin(xi)*sin(xi)));

    ans = a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta;

  } else {

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

    double error;

    gsl_function F;

    double *v = malloc(2*sizeof(double));

    v[0] = xi;
    v[1] = zeta;

    F.function = &f_2_integrand;
    F.params = v;

    gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);

    gsl_integration_workspace_free (w);

    free(v);

  }

  return ans;

}

double g_1func(double x_1,double x_2,double xi,double zeta)
{

  double g_1_integrand(double u,void *params)

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

  double error;

  gsl_function F;

  double *v = malloc(2*sizeof(double));

  v[0] = xi;
  v[1] = zeta;
  
  F.function = &g_1_integrand;
  F.params = v;

  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);

  free(v);

  return ans;

}

double g_2func(double x_1,double x_2,double xi,double zeta)
{

  double g_2_integrand(double u,void *params)

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

  double error;

  gsl_function F;

  double *v = malloc(2*sizeof(double));

  v[0] = xi;
  v[1] = zeta;
  
  F.function = &g_2_integrand;
  F.params = v;

  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);

  free(v);

  return ans;

}

double f_1_integrand(double u,void *params)
// params is just a vector with x[0] = xi, x[1] = zeta.
{
  double *x = *params;

  double f;

  if (fabs(u) < ZERO) {

    f = (2*zeta)*(2*zeta)*u;

  } else {

    f = sin(2*(x[1]*u+x[0]))*sin(2*(x[1]*u+x[0]))/u;

  }

  return f;

}


double f_2_integrand(double u,void *params)
// params is just a vector with x[0] = xi, x[1] = zeta.
{
  double *x = *params;

  double f;

  if (fabs(u) < ZERO) {

    f = x[1]*x[1]*x[1]*x[1]*u*u*u;

  } else {

    f = sin(x[1]*u+x[0])*sin(x[1]*u+x[0])*sin(x[1]*u+x[0])*sin(x[1]*u+x[0])/u;

  }

  return f;

}

double g_1_integrand(double u,void *params)
// params is just a vector with x[0] = xi, x[1] = zeta.
{
  double f;

  f = u/(cos(x[1]*u+x[0])*cos(x[1]*u+x[0]));

  return f;

}

double g_2_integrand(double u,void *params)
// params is just a vector with x[0] = xi, x[1] = zeta.
{
  double f;

  f = u/(cos(x[1]*u+x[0])*cos(x[1]*u+x[0])*cos(x[1]*u+x[0])*cos(x[1]*u+x[0]));

  return f;

}

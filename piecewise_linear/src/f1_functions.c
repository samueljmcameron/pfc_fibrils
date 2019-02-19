#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "headerfile.h"


double f_1func(double x_1,double x_2, double xi,double zeta)
{

  double f_1_integrand(double u,void *params);

  double ans;

  if (fabs(zeta)<SMALL) {

    double a0,a1,a2,a3;

    if (x_1 <= ZERO && xi <= ZERO) {
      
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



double df_1dx_1(double x_1,double xi,double zeta)
{

  double f;

  if (fabs(x_1)<SMALL) {
    
    f = -(2*zeta)*(2*zeta)*x_1;

  } else {

    f = -sin(2*(zeta*x_1+xi))*sin(2*(zeta*x_1+xi))/x_1;

  }

  return f;
}

double df_1dx_2(double x_2,double xi,double zeta)
{
  double df_1dx_1(double x_1,double xi,double zeta);

  return -1*df_1dx_1(x_2,xi,zeta);
}

double df_1dxi(double x_1,double x_2,double xi,double zeta)
{

  double df_1dxi_integrand(double u, void *params);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &df_1dxi_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);
  
  return ans;

}


double df_1dzeta(double x_1,double x_2,double xi,double zeta)
{

  double f;

  if (fabs(zeta)<ZERO) {

    f = 4*(x_2-x_1)*cos(2*xi)*sin(2*xi);

    f += 4*zeta*(x_2*x_2-x_1*x_1)*(cos(2*xi)*cos(2*xi)-sin(2*xi)*sin(2*xi));

  } else {

    f =  1/(zeta)*(sin(2*(zeta*x_2+xi))*sin(2*(zeta*x_2+xi))
		   -sin(2*(zeta*x_1+xi))*sin(2*(zeta*x_1+xi)));

  }

  return f;

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



double df_1dxi_integrand(double u, void *params)
{
  double *x = params;

  x[0] = xi;
  x[1] = zeta;

  double f;

  if (fabs(u) < ZERO) {

    f = 4*(2*zeta);

  } else {

    f = 4*sin(2*(zeta*u+xi))*cos(2*(zeta*u+xi))/u;

  }

  return f;

}


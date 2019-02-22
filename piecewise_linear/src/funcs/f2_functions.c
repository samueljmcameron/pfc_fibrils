#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "../headerfile.h"

double f_2func(double x_1,double x_2, double xi,double zeta)
{

  double f_2_integrand(double u,void *params);

  double ans;

  if (fabs(zeta)<SMALL) {

    double a0,a1,a2,a3;

    if ((x_1 < ZERO || x_2 < ZERO) && fabs(xi) < ZERO) {
      
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

    gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,w,&ans,&error);

    gsl_integration_workspace_free (w);

    free(v);

  }



  return ans;

}


double df_2dx_1(double x_1,double xi,double zeta)
{

  double f;

  if (fabs(x_1)<SMALL) {
    
    f = -(zeta)*(zeta)*x_1*x_1*x_1;

  } else {

    f = -sin(zeta*x_1+xi)*sin(zeta*x_1+xi)*sin(zeta*x_1+xi)*sin(zeta*x_1+xi)/x_1;

  }

  return f;
}

double df_2dx_2(double x_2,double xi,double zeta)
{
  double df_2dx_1(double x_1,double xi,double zeta);

  return -1*df_2dx_1(x_2,xi,zeta);
}


double df_2dxi(double x_1,double x_2,double xi,double zeta)
{

  double df_2dxi_integrand(double u, void *params);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

  double ans;
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &df_2dxi_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,w,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);

  return ans;

}

double df_2dzeta(double x_1,double x_2,double xi,double zeta)
{

  double f;

  if (fabs(zeta)<SMALL) {

    f = 4*(x_2-x_1)*sin(xi)*sin(xi)*sin(xi)*cos(xi);

    f += 2*zeta*(x_2*x_2-x_1*x_1)*sin(xi)*sin(xi)*(3*cos(xi)*cos(xi)-sin(xi)*sin(xi));          
    
  } else {

    f = 1/(zeta)*(sin(zeta*x_2+xi)*sin(zeta*x_2+xi)*sin(zeta*x_2+xi)*sin(zeta*x_2+xi)
		  -sin(zeta*x_1+xi)*sin(zeta*x_1+xi)*sin(zeta*x_1+xi)*sin(zeta*x_1+xi));

  }

  return f;
}


double f_2_integrand(double u,void *params)
// params is just a vector with x[0] = xi, x[1] = zeta.
{
  double *x = params;

  double xi = x[0];
  double zeta = x[1];

  double f;

  if (fabs(u) < ZERO) {

    f = zeta*zeta*zeta*zeta*u*u*u;

  } else {

    f = sin(zeta*u+xi)*sin(zeta*u+xi)*sin(zeta*u+xi)*sin(zeta*u+xi)/u;

  }

  return f;

}



double df_2dxi_integrand(double u, void *params)
{

  double *x = params;

  double xi = x[0];
  double zeta = x[1];


  double f;

  if (fabs(u) < ZERO) {

    f = 4*zeta*zeta*zeta*u*u;

  } else {

    f = 4*sin(zeta*u+xi)*sin(zeta*u+xi)*sin(zeta*u+xi)*cos(zeta*u+xi)/u;

  }

  return f;

}


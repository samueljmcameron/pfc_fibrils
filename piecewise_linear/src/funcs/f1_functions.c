#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "../headerfile.h"


double f_1func(double x_1,double x_2, double xi,double zeta)
{

  double f_1_integrand(double u,void *params);

  double ans;

  if (fabs(zeta)<SMALL) {

    double x_2t2 = x_2*x_2;

    double x_2t3 = x_2*x_2*x_2;

    double x_2t4 = x_2*x_2*x_2*x_2;

    double x_2t5 = x_2*x_2*x_2*x_2*x_2;

    double x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2;

    double x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2;

    double x_1t2 = x_1*x_1;

    double x_1t3 = x_1*x_1*x_1;

    double x_1t4 = x_1*x_1*x_1*x_1;

    double x_1t5 = x_1*x_1*x_1*x_1*x_1;

    double x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1;

    double x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1;


    double a0,a1,a2,a3,a4,a5;

    if (x_1 <= ZERO && xi <= ZERO) {
      
      a0 = 0.0;

    } else {
      
      a0 = sin(2*xi)*sin(2*xi)*log(x_2/x_1);

    }
    /*
    a1 = 4*(x_2-x_1)*cos(2*xi)*sin(2*xi);

    a2 = 2*(x_2*x_2-x_1*x_1)*(cos(2*xi)*cos(2*xi)-sin(2*xi)*sin(2*xi));

    a3 = -32.0/9.0*(x_2*x_2*x_2-x_1*x_1*x_1)*sin(2*xi)*cos(2*xi);

    */

    double cos2xi_t2 = cos(2*xi)*cos(2*xi);

    double sin2xi_t2 = sin(2*xi)*sin(2*xi);

    a1 = 4*sin(2*xi)*cos(2*xi)*(x_2-x_1);

    a2 = +(-2*cos2xi_t2*x_1t2+2*cos2xi_t2*x_2t2+2*sin2xi_t2*x_1t2
	   -2*sin2xi_t2*x_2t2);

    a3 =-(32/9.0)*sin(2*xi)*cos(2*xi)*(-x_1t3+x_2t3);

    a4 = +(4*cos2xi_t2*x_1t4*(1/3.0)-4*sin2xi_t2*x_1t4*(1/3.0)
	   -4*cos2xi_t2*x_2t4*(1/3.0)+4*sin2xi_t2*x_2t4*(1/3.0));

    a5 = +(128/75.0)*sin(2*xi)*cos(2*xi)*(-x_1t5+x_2t5);

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);

  } else {

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);

    double error;

    gsl_function F;

    double *v = malloc(2*sizeof(double));

    v[0] = xi;
    v[1] = zeta;

    F.function = &f_1_integrand;
    F.params = v;

    gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,w,&ans,&error);
    
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

  double ans;
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &df_1dxi_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,w,&ans,&error);
  
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
  
  double *x = params;
  

  double xi=x[0];
  double zeta=x[1];

  double f;

  if (fabs(u) < ZERO) {

    f = (2*zeta)*(2*zeta)*u;

  } else {

    f = sin(2*(zeta*u+xi))*sin(2*(zeta*u+xi))/u;

  }

  return f;

}



double df_1dxi_integrand(double u, void *params)
{
  double *x = params;

  double xi=x[0];
  double zeta=x[1];
  
  double f;

  if (fabs(u) < ZERO) {

    f = 4*(2*zeta);

  } else {

    f = 4*sin(2*(zeta*u+xi))*cos(2*(zeta*u+xi))/u;

  }

  return f;

}


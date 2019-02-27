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

    if ((x_1 < ZERO || x_2 < ZERO) && fabs(xi) < ZERO) {
      
      a0 = 0.0;

    } else {
      
      a0 = sin(xi)*sin(xi)*sin(xi)*sin(xi)*log(x_2/x_1);

    }

    double sinxi_t2 = sin(xi)*sin(xi);

    double sinxi_t3 = sin(xi)*sinxi_t2;

    double sinxi_t4 = sinxi_t3*sin(xi);

    double cosxi_t2 = cos(xi)*cos(xi);

    double cosxi_t3 = cos(xi)*cosxi_t2;

    double cosxi_t4 = cosxi_t3*cos(xi);


    /*
    a1 = 4*(x_2-x_1)*cos(xi)*sin(xi)*sin(xi)*sin(xi);

    a2 = (x_2*x_2-x_1*x_1)*sin(xi)*sin(xi)*(3*cos(xi)*cos(xi)-sin(xi)*sin(xi));

    a3 = (4.0/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*sin(xi)*cos(xi)
	  *(cos(xi)*cos(xi)-5*sin(xi)*sin(xi)));
	  */


    a1 = +4*sinxi_t3*cos(xi)*(x_2-x_1);

    a2 = +(sinxi_t4*x_1t2-sinxi_t4*x_2t2-3*sinxi_t2*cosxi_t2*x_1t2
	   +3*sinxi_t2*cosxi_t2*x_2t2);

    a3 = +((20*sinxi_t3*cos(xi)*x_1t3)/9.0-(4*sin(xi)*cosxi_t3*x_1t3)/3.0
	   -(20*sinxi_t3*cos(xi)*x_2t3)/9.0+(4*sin(xi)*cosxi_t3*x_2t3)/3.0);

    a4 = +(-(5*sinxi_t4*x_1t4)/12.0+2*sinxi_t2*cosxi_t2*x_1t4
	   -(cosxi_t4*x_1t4)/4.0+(5*sinxi_t4*x_2t4)/12.0
	   -2*sinxi_t2*cosxi_t2*x_2t4+(cosxi_t4*x_2t4)/4.0);

    a5 = +(-(68*sinxi_t3*cos(xi)*x_1t5)/75.0+(4*sin(xi)*cosxi_t3*x_1t5)/5.0
	   +(68*sinxi_t3*cos(xi)*x_2t5)/75.0-(4*sin(xi)*cosxi_t3*x_2t5)/5.0);

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);

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


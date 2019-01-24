#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "headerfile.h"

double dudx_1(double x_1,double zeta)
{

  return -2*(1-zeta)*(1-zeta)*x_1;

}

double dudx_2(double x_2,double zeta)
{

  return 2*(1-zeta)*(1-zeta)*x_2;

}

double dudxi(void)
{

  return 0;

}

double dudzeta(double x_1,double x_2,double zeta)
{

  return -2*zeta*(1-zeta)*(x_2*x_2-x_1*x_1);

}


double dvdx_1(double x_1,double xi,double zeta)
{

  return 2*(1-zeta)*sin(2*(zeta*x_1+xi));

}

double dvdx_2(double x_2,double xi,double zeta)
{

  return -2*(1-zeta)*sin(2*(zeta*x_2+xi));

}

double dvdxi(double x_1,double x_2,double xi, double zeta)
{

  double f;

  if (fabs(zeta)<SMALL) {

    f = -4*cos(2*xi)*(x_2-x_1);
    f += zeta*(4*cos(2*xi)*(x_2-x_1)+4*sin(2*xi)*(x_2*x_2-x_1*x_1));

  } else {

    f = -2*(1-zeta)/zeta*(sin(2*(zeta*x_2+xi))-sin(2*(zeta*x_1+xi)));

  }

  return f;

}

double dvdzeta(double x_1,double x_2,double xi,double zeta)
{
  double f;

  if (fabs(zeta)<SMALL) {

    f = 2*sin(2*xi)*(x_2-x_1)-2*cos(2*xi)*(x_2*x_2-x_1*x_1);
    f += 4*zeta*(cos(2*xi)*(x_2*x_2-x_1*x_1)
		 +2.0/3.0*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1));

  } else {

    f = (1/zeta-1)*(2*x_1*sin(2*(zeta*x_1+xi))-2*x_2*sin(2*(zeta*x_2+xi)));

    f += -1/(zeta*zeta)*(np.cos(2*(zeta*x_2+xi))-np.cos(2*(zeta*x_1+xi)));

  }

  return f;

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

double df_2dxi(double x_1,double x_2,double xi,double zeta)
{

  double df_2dxi_integrand(double u, void *params);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &df_2dxi_integrand;
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

double dg_1dx_1(double x_1,double xi,double zeta)
{
  double f;

  f = -x_1/(cos(zeta*x_1+xi)*cos(zeta*x_1+xi));

  return f;

}

double dg_1dx_2(double x_2,double xi,double zeta)
{
 
  return -dg_1dx_1(x_2,xi,zeta);
  
}

double dg_1dxi(double x_1,double x_2,double xi,double zeta)
{

  double dg_1dxi_integrand(double u,void *params);


  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &dg_1dxi_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);

  return ans;


}

double dg_2dxi(double x_1,double x_2,double xi,double zeta)
{

  double dg_2dxi_integrand(double u,void *params);


  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &dg_2dxi_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);

  return ans;


}

double dg_1dzeta(double x_1,double x_2,double xi,double zeta)
{

  double dg_1dzeta_integrand(double u,void *params);


  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &dg_1dzeta_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);

  return ans;


}

double dg_2dzeta(double x_1,double x_2,double xi,double zeta)
{

  double dg_2dzeta_integrand(double u,void *params);


  gsl_integration_workspace *w = gsl_integration_workspace_alloc(SIZE_INTEGRATION);
  
  double error;
  
  gsl_function F;
  
  double *v = malloc(2*sizeof(double));
  
  v[0] = xi;
  v[1] = zeta;
  
  F.function = &dg_2dzeta_integrand;
  F.params = v;
  
  gsl_integration_qags(&F,x_1,x_2,ABSERR,RELERR,SIZE_INTEGRATION,&ans,&error);
  
  gsl_integration_workspace_free (w);
  
  free(v);

  return ans;


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

double df_2dxi_integrand(double u, void *params)
{

  double *x = params;

  x[0] = xi;
  x[1] = zeta;

  double f;

  if (fabs(u) < ZERO) {

    f = 4*zeta*zeta*zeta*u*u;

  } else {

    f = 4*sin(zeta*u+xi)*sin(zeta*u+xi)*sin(zeta*u+xi)*cos(zeta*u+xi)/u;

  }

  return f;

}


double dg_1dxi_integrand(double u,void *params)
{

  double *x = params;

  x[0] = xi;
  x[1] = zeta;

  double f;

  f = 2*u*sin(zeta*u+xi)/(cos(zeta*u+xi)*cos(zeta*u+xi)*cos(zeta*u+xi));

  return f;

}

double dg_2dxi_integrand(double u,void *params)
{

  double *x = params;

  x[0] = xi;

  x[1] = zeta;

  double f;

  f = 4*u*sin(zeta*u+xi)/(cos(zeta*u+xi)*cos(zeta*u+xi)*cos(zeta*u+xi)
			  *cos(zeta*u+xi)*cos(zeta*u+xi));

  return f;

}

double dg_1dzeta_integrand(double u,void *params)
{

  double *x = params;

  x[0] = xi;
  x[1] = zeta;

  double f;

  f = 2*u*u*sin(zeta*u+xi)/(cos(zeta*u+xi)*cos(zeta*u+xi)*cos(zeta*u+xi));

  return f;

}

double dg_2dzeta_integrand(double u,void *params)
{

  double *x = params;

  x[0] = xi;

  x[1] = zeta;

  double f;

  f = 4*u*u*sin(zeta*u+xi)/(cos(zeta*u+xi)*cos(zeta*u+xi)*cos(zeta*u+xi)
			    *cos(zeta*u+xi)*cos(zeta*u+xi));

  return f;

}

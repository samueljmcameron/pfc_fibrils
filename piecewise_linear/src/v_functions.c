#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "headerfile.h"

double vfunc(double x_1,double x_2,double xi,double zeta)
{
 
  double ans;
  
  if (fabs(zeta)< SMALL) {

    double a0,a1,a2,a3;

    a0 = -2*sin(2*xi)*(x_2-x_1);

    a1 = 2*(sin(2*xi)*(x_2-x_1)-cos(2*xi)*(x_2*x_2-x_1*x_1));

    a2 = 2.0/3.0*(3*cos(2*xi)*(x_2*x_2-x_1*x_1)+2*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1));

    a3 = 2.0/3.0*(-2*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1)
		  +cos(2*xi)*(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1));
    
    ans = a0+a1*zeta+a2*zeta*zeta+a3*zeta*zeta*zeta;
    
  } else {

    ans = cos(2*(zeta*x_2+xi))-cos(2*(zeta*x_1+xi));

    ans *= (1-zeta)/zeta;

  }

  return ans;

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



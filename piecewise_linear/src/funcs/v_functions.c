#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"

double vfunc(double x_1,double x_2,double xi,double zeta)
{
 
  double ans;
  
  if (fabs(zeta)< SMALL) {

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

    //a0 = -2*sin(2*xi)*(x_2-x_1);

    //a1 = 2*(sin(2*xi)*(x_2-x_1)-cos(2*xi)*(x_2*x_2-x_1*x_1));

    //a2 = 2.0/3.0*(3*cos(2*xi)*(x_2*x_2-x_1*x_1)+2*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1));

    //a3 = 2.0/3.0*(-2*sin(2*xi)*(x_2*x_2*x_2-x_1*x_1*x_1)
    //	  +cos(2*xi)*(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1));

    a0 = -2*sin(2*xi)*x_2+2*sin(2*xi)*x_1;

    a1 = +(-2*cos(2*xi)*x_2t2+2*cos(2*xi)*x_1t2+2*sin(2*xi)*x_2-2*sin(2*xi)*x_1);

    a2 = +(4*sin(2*xi)*x_2t3*(1/3.0)-4*sin(2*xi)*x_1t3*(1/3.0)
	   +2*cos(2*xi)*x_2t2-2*cos(2*xi)*x_1t2);

    a3 = +(2*cos(2*xi)*x_2t4*(1/3.0)-2*cos(2*xi)*x_1t4*(1/3.0)
	   -4*sin(2*xi)*x_2t3*(1/3.0)+4*sin(2*xi)*x_1t3*(1/3.0));

    a4 = +(-4*sin(2*xi)*x_2t5*(1/15.0)+4*sin(2*xi)*x_1t5*(1/15.0)
	   -2*cos(2*xi)*x_2t4*(1/3.0)+2*cos(2*xi)*x_1t4*(1/3.0));

    a5 = +(-4*cos(2*xi)*x_2t6*(1/45.0)+4*cos(2*xi)*x_1t6*(1/45.0)
	   +4*sin(2*xi)*x_2t5*(1/15.0)-4*sin(2*xi)*x_1t5*(1/15.0));

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);
    
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

    a0 = -4*cos(2*xi)*x_2+4*cos(2*xi)*x_1;

    a1 = +(4*sin(2*xi)*x_2t2-4*sin(2*xi)*x_1t2+4*cos(2*xi)*x_2-4*cos(2*xi)*x_1);

    a2 = +(8*cos(2*xi)*x_2t3*(1/3.0)-8*cos(2*xi)*x_1t3*(1/3.0)
	   -4*sin(2*xi)*x_2t2+4*sin(2*xi)*x_1t2);

    a3 = +(-4*sin(2*xi)*x_2t4*(1/3.0)+4*sin(2*xi)*x_1t4*(1/3.0)
	   -8*cos(2*xi)*x_2t3*(1/3.0)+8*cos(2*xi)*x_1t3*(1/3.0));

    a4 = +(-8*cos(2*xi)*x_2t5*(1/15.0)+8*cos(2*xi)*x_1t5*(1/15.0)
	   +4*sin(2*xi)*x_2t4*(1/3.0)-4*sin(2*xi)*x_1t4*(1/3.0));

    a5 = +(8*sin(2*xi)*x_2t6*(1/45.0)-8*sin(2*xi)*x_1t6*(1/45.0)
	   +8*cos(2*xi)*x_2t5*(1/15.0)-8*cos(2*xi)*x_1t5*(1/15.0));

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);


  } else {

    ans = -2*(1-zeta)/zeta*(sin(2*(zeta*x_2+xi))-sin(2*(zeta*x_1+xi)));

  }

  return ans;

}

double dvdzeta(double x_1,double x_2,double xi,double zeta)
{
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

    a0 = -2*cos(2*xi)*x_2t2+2*cos(2*xi)*x_1t2+2*sin(2*xi)*x_2-2*sin(2*xi)*x_1;

    a1 = +(8*sin(2*xi)*x_2t3*(1/3.0)-8*sin(2*xi)*x_1t3*(1/3.0)
	   +4*cos(2*xi)*x_2t2-4*cos(2*xi)*x_1t2);
    
    a2 = +(2*cos(2*xi)*x_2t4-2*cos(2*xi)*x_1t4
	   -4*sin(2*xi)*x_2t3+4*sin(2*xi)*x_1t3);
    
    a3 = +(-16*sin(2*xi)*x_2t5*(1/15.0)+16*sin(2*xi)*x_1t5*(1/15.0)
	   -8*cos(2*xi)*x_2t4*(1/3.0)+8*cos(2*xi)*x_1t4*(1/3.0));
    
    a4 = +(-4*cos(2*xi)*x_2t6*(1/9.0)+4*cos(2*xi)*x_1t6*(1/9.0)
	   +4*sin(2*xi)*x_2t5*(1/3.0)-4*sin(2*xi)*x_1t5*(1/3.0));
    
    a5 = +(16*sin(2*xi)*x_2t7*(1/105.0)-16*sin(2*xi)*x_1t7*(1/105.0)
	   +8*cos(2*xi)*x_2t6*(1/15.0)-8*cos(2*xi)*x_1t6*(1/15.0));

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);
    
    
  } else {

    ans = (1/zeta-1)*(2*x_1*sin(2*(zeta*x_1+xi))-2*x_2*sin(2*(zeta*x_2+xi)));

    ans += -1/(zeta*zeta)*(cos(2*(zeta*x_2+xi))-cos(2*(zeta*x_1+xi)));

  }

  return ans;

}



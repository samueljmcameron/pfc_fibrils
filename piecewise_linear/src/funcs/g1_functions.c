#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"

double g_1func(double x_1,double x_2,double xi,double zeta)
{

  double ans;
  
  if (fabs(zeta)< SMALL) {

    double a0,a1,a2,a3;

    a0 = (x_2*x_2-x_1*x_1)*cos(xi)*cos(xi)/2.0;
      
    a1 = -2*(x_2*x_2*x_2-x_1*x_1*x_1)*cos(xi)*sin(xi)/3.0;

    a2 = -(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1)*(cos(xi)*cos(xi)-sin(xi)*sin(xi))/4.0;

    a3 = 4*(x_2*x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1*x_1)*cos(xi)*sin(xi)/15.0;

    ans = a0+a1*zeta+a2*zeta*zeta+a3*zeta*zeta*zeta;

  } else {

    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);
    
    ans = 2*zeta*(x_2*sin(s2)-x_1*sin(s1))+cos(s2)-cos(s1);
    ans += 2*zeta*zeta*(x_2*x_2-x_1*x_1);
      
    ans /= 8*zeta*zeta;

  }


  return ans;

}


double dg_1dx_1(double x_1,double xi,double zeta)
{
  double f;

  f = -x_1*(cos(zeta*x_1+xi)*cos(zeta*x_1+xi));

  return f;

}

double dg_1dx_2(double x_2,double xi,double zeta)
{
 
  return -dg_1dx_1(x_2,xi,zeta);
  
}

double dg_1dxi(double x_1,double x_2,double xi,double zeta)
{

  double ans;
  
  if (fabs(zeta)<SMALL) {

    double a0,a1,a2,a3;

    a0 = -(x_2*x_2-x_1*x_1)*cos(xi)*sin(xi);

    a1 = -2/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*(cos(xi)*cos(xi)-sin(xi)*sin(xi));

    a2 = 0.5*(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1)*cos(xi)*sin(xi);

    a3 = 1./15.*((x_2*x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1*x_1)
		 *(20*cos(xi)*cos(xi)*cos(xi)*sin(xi)
		   +12*sin(xi)*sin(xi)*(3*cos(xi)*cos(xi)-sin(xi)*sin(xi))));

    ans = a0 + a1*zeta;// + a2*zeta*zeta + a3*zeta*zeta*zeta;

  } else {

    double p1 = zeta*x_1+xi;
    double p2 = zeta*x_2+xi;
        
    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);

    ans = sin(s2)-sin(s1)-(s2*cos(s2)-s1*cos(s1));
    ans += 4*xi*(cos(p2)*cos(p2)-cos(p1)*cos(p1));
      
    ans /= 4*zeta*zeta;

    ans *= -1;

  }

  return ans;

}

double dg_1dzeta(double x_1,double x_2,double xi,double zeta)
{

  double ans;


  if (fabs(zeta)<SMALL) {

    double a0,a1,a2;

    a0 = -2/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*cos(xi)*sin(xi);
            
    a1 = -1/2.0*(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1)*(cos(xi)*cos(xi)-sin(xi)*sin(xi));

    a2 = 4*(x_2*x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1*x_1)*cos(xi)*sin(xi)/5.0;

    ans = a0 + a1*zeta;//+a2*zeta*zeta;

  } else {

    double p1 = zeta*x_1+xi;
    double p2 = zeta*x_2+xi;
        
    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);

    double q1 = -2*x_1*x_1*zeta*zeta+2*xi*xi+1;
    double q2 = -2*x_2*x_2*zeta*zeta+2*xi*xi+1;
            
    ans = 2*zeta*(x_2*sin(s2)-x_1*sin(s1));
    ans += q2*cos(s2)-q1*cos(s1);
    ans += -4*xi*xi*(cos(p2)*cos(p2)-cos(p1)*cos(p1));

    ans /= 4*zeta*zeta*zeta;

    ans *= -1;

  }

  return ans;

}


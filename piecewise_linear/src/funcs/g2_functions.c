#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"

double g_2func(double x_1,double x_2,double xi,double zeta)
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
    
    a0 = (x_2*x_2-x_1*x_1)*cos(xi)*cos(xi)*cos(xi)*cos(xi)/2.0;

    a1 = -4*(x_2*x_2*x_2-x_1*x_1*x_1)*cos(xi)*cos(xi)*cos(xi)*sin(xi)/3.0;
    
    a2 = -((x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1)*cos(xi)*cos(xi)*
	   (cos(xi)*cos(xi)-3*sin(xi)*sin(xi))/2.0);
    
    a3 = -((x_2*x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1*x_1)*cos(xi)
	   *(5*cos(xi)*cos(xi)*cos(xi)-12*sin(xi)*sin(xi)*sin(xi))/15.0);
    
    
    a0 = ((1/16.0)*cos(4*xi)*x_2t2+(1/4.0)*cos(2*xi)*x_2t2-(1/16.0)*cos(4*xi)*x_1t2
	  -(1/4.0)*cos(2*xi)*x_1t2-3*x_1t2*(1/16.0)+3*x_2t2*(1/16.0));

    a1 = (-(1/6.0)*sin(4*xi)*x_2t3-(1/3.0)*sin(2*xi)*x_2t3
	  +(1/6.0)*sin(4*xi)*x_1t3+(1/3.0)*sin(2*xi)*x_1t3);

    a2 = (-(1/4.0)*cos(4*xi)*x_2t4-(1/4.0)*cos(2*xi)*x_2t4
	  +(1/4.0)*cos(4*xi)*x_1t4+(1/4.0)*cos(2*xi)*x_1t4);
            
    a3 = (4*sin(4*xi)*x_2t5*(1/15.0)+2*sin(2*xi)*x_2t5*(1/15.0)
	  -4*sin(4*xi)*x_1t5*(1/15.0)-2*sin(2*xi)*x_1t5*(1/15.0));

    a4 = (2*cos(4*xi)*x_2t6*(1/9.0)+(1/18.0)*x_2t6*cos(2*xi)
	  -2*cos(4*xi)*x_1t6*(1/9.0)-(1/18.0)*x_1t6*cos(2*xi));

    a5 = (-16*sin(4*xi)*x_2t7*(1/105.0)-2*sin(2*xi)*x_2t7*(1/105.0)
    	  +16*sin(4*xi)*x_1t7*(1/105.0)+2*sin(2*xi)*x_1t7*(1/105.0));
    

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);


  } else {

    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);

    
    double t1 = 4*(zeta*x_1+xi);
    double t2 = 4*(zeta*x_2+xi);

    ans = 4*zeta*(x_2*sin(t2)-x_1*sin(t1))+cos(t2)-cos(t1);
    ans += 32*zeta*(x_2*sin(s2)-x_1*sin(s1));
    ans += 16*(cos(s2)-cos(s1))+24*zeta*zeta*(x_2*x_2-x_1*x_1);

    ans /= 128*zeta*zeta;

  }

  return ans;

}


double dg_2dx_1(double x_1,double xi,double zeta)
{

  double f;

  f = -x_1*cos(zeta*x_1+xi)*cos(zeta*x_1+xi)*cos(zeta*x_1+xi)*cos(zeta*x_1+xi);

  return f;

}

double dg_2dx_2(double x_2,double xi,double zeta)
{
  
  double dg_2dx_1(double x_1,double xi,double zeta);
 
  return -dg_2dx_1(x_2,xi,zeta);
  
}


double dg_2dxi(double x_1,double x_2,double xi,double zeta)
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

    a0 = -2*(x_2*x_2-x_1*x_1)*cos(xi)*cos(xi)*cos(xi)*sin(xi);

    a1 = -4/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*cos(xi)*cos(xi)*(cos(xi)*cos(xi)-3*sin(xi)*sin(xi));

    a0 = (-(1/4.0)*sin(4*xi)*x_2t2-(1/2.0)*sin(2*xi)*x_2t2+(1/4.0)*sin(4*xi)*x_1t2
	  +(1/2.0)*sin(2*xi)*x_1t2);

    a1 = (-2*cos(4*xi)*x_2t3*(1/3.0)-2*cos(2*xi)*x_2t3*(1/3.0)
	  +2*cos(4*xi)*x_1t3*(1/3.0)+2*cos(2*xi)*x_1t3*(1/3.0));

    a2 = (sin(4*xi)*x_2t4+(1/2.0)*sin(2*xi)*x_2t4-sin(4*xi)*x_1t4
	  -(1/2.0)*sin(2*xi)*x_1t4);

    a3 = +(16*cos(4*xi)*x_2t5*(1/15.0)+4*cos(2*xi)*x_2t5*(1/15.0)
	   -16*cos(4*xi)*x_1t5*(1/15.0)-4*cos(2*xi)*x_1t5*(1/15.0))*0;

    a4 = +(-8*sin(4*xi)*x_2t6*(1/9.0)-(1/9.0)*sin(2*xi)*x_2t6
	   +8*sin(4*xi)*x_1t6*(1/9.0)+(1/9.0)*sin(2*xi)*x_1t6)*0;

    a5 = +(-64*cos(4*xi)*x_2t7*(1/105.0)-4*cos(2*xi)*x_2t7*(1/105.0)
	   +64*cos(4*xi)*x_1t7*(1/105.0)+4*cos(2*xi)*x_1t7*(1/105.0))*0;


    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta + a4*zeta*zeta*zeta*zeta
	   + a5*zeta*zeta*zeta*zeta*zeta);

  } else {

    double p1 = zeta*x_1+xi;
    double p2 = zeta*x_2+xi;
        
    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);

    double t1 = 4*(zeta*x_1+xi);
    double t2 = 4*(zeta*x_2+xi);

    ans = sin(t2)-sin(t1)-(t2*cos(t2)-t1*cos(t1));
    ans += 8*(sin(s2)-sin(s1))-8*(s2*cos(s2)-s1*cos(s1));
    ans += 32*xi*(cos(p2)*cos(p2)*cos(p2)*cos(p2)-cos(p1)*cos(p1)*cos(p1)*cos(p1));

    ans /= 32*zeta*zeta;

    ans *= -1;

  }

  return ans;


}


double dg_2dzeta(double x_1,double x_2,double xi,double zeta)
{

  double ans;

  if (fabs(zeta)<SMALL) {

    double x_2t2 = x_2*x_2;

    double x_2t3 = x_2*x_2*x_2;

    double x_2t4 = x_2*x_2*x_2*x_2;

    double x_2t5 = x_2*x_2*x_2*x_2*x_2;

    double x_2t6 = x_2*x_2*x_2*x_2*x_2*x_2;

    double x_2t7 = x_2*x_2*x_2*x_2*x_2*x_2*x_2;

    double x_2t8 = x_2*x_2*x_2*x_2*x_2*x_2*x_2*x_2;

    double x_1t2 = x_1*x_1;

    double x_1t3 = x_1*x_1*x_1;

    double x_1t4 = x_1*x_1*x_1*x_1;

    double x_1t5 = x_1*x_1*x_1*x_1*x_1;

    double x_1t6 = x_1*x_1*x_1*x_1*x_1*x_1;

    double x_1t7 = x_1*x_1*x_1*x_1*x_1*x_1*x_1;

    double x_1t8 = x_1*x_1*x_1*x_1*x_1*x_1*x_1*x_1;


    double a0,a1,a2,a3,a4,a5;

    a0 = -4/3.0*(x_2*x_2*x_2-x_1*x_1*x_1)*cos(xi)*cos(xi)*cos(xi)*sin(xi);
            
    a1 = -(x_2*x_2*x_2*x_2-x_1*x_1*x_1*x_1)*cos(xi)*cos(xi)*(cos(xi)*cos(xi)-sin(xi)*sin(xi));


    a0 = (-(1/6.0)*sin(4*xi)*x_2t3-(1/3.0)*sin(2*xi)*x_2t3
	  +(1/6.0)*sin(4*xi)*x_1t3+(1/3.0)*sin(2*xi)*x_1t3);

    a1 = (-(1/2.0)*cos(4*xi)*x_2t4-(1/2.0)*cos(2*xi)*x_2t4+(1/2.0)*cos(4*xi)*x_1t4
	  +(1/2.0)*cos(2*xi)*x_1t4);

    a2 = (4*sin(4*xi)*x_2t5*(1/5.0)+2*sin(2*xi)*x_2t5*(1/5.0)
	  -4*sin(4*xi)*x_1t5*(1/5.0)-2*sin(2*xi)*x_1t5*(1/5.0));

    a3 = (8*cos(4*xi)*x_2t6*(1/9.0)+2*x_2t6*cos(2*xi)*(1/9.0)
	  -8*cos(4*xi)*x_1t6*(1/9.0)-2*x_1t6*cos(2*xi)*(1/9.0))*0;

    a4 = (-16*sin(4*xi)*x_2t7*(1/21.0)-2*sin(2*xi)*x_2t7*(1/21.0)
	  +16*sin(4*xi)*x_1t7*(1/21.0)+2*sin(2*xi)*x_1t7*(1/21.0))*0;

    a5 = (-8*cos(4*xi)*x_2t8*(1/15.0)-(1/30.0)*x_2t8*cos(2*xi)
	  +8*cos(4*xi)*x_1t8*(1/15.0)+(1/30.0)*x_1t8*cos(2*xi))*0;

    ans = (a0 + a1*zeta + a2*zeta*zeta + a3*zeta*zeta*zeta
	   + a4*zeta*zeta*zeta*zeta + a5*zeta*zeta*zeta*zeta*zeta);

  } else {

    double p1 = zeta*x_1+xi;
    double p2 = zeta*x_2+xi;
        
    double s1 = 2*(zeta*x_1+xi);
    double s2 = 2*(zeta*x_2+xi);

    double q1 = -2*x_1*x_1*zeta*zeta+2*xi*xi+1;
    double q2 = -2*x_2*x_2*zeta*zeta+2*xi*xi+1;

    double t1 = 4*(zeta*x_1+xi);
    double t2 = 4*(zeta*x_2+xi);

    double r1 = -8*x_1*x_1*zeta*zeta+8*xi*xi+1;
    double r2 = -8*x_2*x_2*zeta*zeta+8*xi*xi+1;

    ans = 4*zeta*(x_2*sin(t2)-x_1*sin(t1))+r2*cos(t2)-r1*cos(t1);
    ans += 32*zeta*(x_2*sin(s2)-x_1*sin(s1));
    ans += 16*(q2*cos(s2)-q1*cos(s1));
    ans += -64*xi*xi*(cos(p2)*cos(p2)*cos(p2)*cos(p2)-cos(p1)*cos(p1)*cos(p1)*cos(p1));

    ans /= 64*zeta*zeta*zeta;

    ans *= -1;
    
  }


  return ans;


}



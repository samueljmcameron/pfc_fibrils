#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"

double ufunc(double x_1,double x_2,double zeta)
{
  return (1-zeta)*(1-zeta)*(x_2*x_2-x_1*x_1);
}

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



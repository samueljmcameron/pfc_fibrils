#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

double trapzd(double *x, double *y, double hmin, int xlength, int n)
/*This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input as a pointer to the function to be integrated between limits a and b, also input. When called with n=1, the routine returns the crudest estimate of integral_a^b(f(x)dx). Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points. */
{
  double tnm,sum,del;
  static double s;
  int it,j,index,spacing;

  if (n == 1) {
    return (s=0.5*(x[xlength]-x[1])*(y[1]+y[xlength]));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it; // equals 1 if n = 2
    del=(x[xlength]-x[1])/tnm; //This is the spacing of the points to be added.
    if (del < 2*hmin) nrerror("Stepsize of integration smaller than discrete data spacing.");
    index = 0+round(del/(2*hmin))+1;
    spacing = 2*(index-1);
    for (sum=0.0,j=1;j<=it;j++,index += spacing) sum += y[index];
    s=0.5*(s+(x[xlength]-x[1])*sum/tnm); //This replaces s by its refined value.
    return s;
  }
}

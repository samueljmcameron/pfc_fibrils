#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "integrate.h"

double trapzd(double *x, double *y, double hmin, int xlength, int n);

int main(void)
{
  int xlength = 2*2*2*2*2*2*2+1;
  double x[xlength];
  double y[xlength];
  double s;
  int j;
  int nmax = round(log((xlength-1))/log(2));
  double hmin;
  for (j = 0; j < xlength; j++) {
    x[j] = j*(M_PI)/(xlength-1);
    y[j] = sin(x[j]);
  }
  hmin = x[2]-x[1];
  nmax = round(log((xlength-1))/log(2));
  for (j=1; j <= nmax + 2; j++) s = trapzd(x-1,y-1,hmin,xlength,j);
  printf("%lf\n",s);
}

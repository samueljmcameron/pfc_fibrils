#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"


#define EPS 1.0e-10
#define K 5

// sample driver integrating sin(x) from 0 to pi
/*
double qromb(double *x,double *y, double hmin, int xlength, int JMAX);

int main(void)
{
  int xlength = 2*2*2*2*2*2+1;
  double x[xlength];
  double y[xlength];
  printf("%d\n",xlength);
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
  s = qromb(x-1,y-1,hmin,xlength,nmax);
  printf("%lf\n",s);
}
*/
/*Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
 JMAX limits the total number of steps; K is the number of points used in the extrapolation.*/

double qromb(double *x,double *y, int xlength,bool *failure)
/*Returns the integral of the function func from a to b. Integration is performed by Romberg's 
method of order 2K, where, e.g., K=2 is Simpsons rule. */
{
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double *x, double *y, double hmin, int xlength, int n);
  void nrerror(char error_text[]);
  double ss=1e50,dss,lastss=1e50;
  int JMAX = round(log((xlength-1))/log(2)); // cannot have more steps than there are data points;
  double s[JMAX+2],h[JMAX+3]; //These store the successive trapezoidal approxi-
  int j;                      //mations and their relative stepsizes.
  double hmin = x[2]-x[1];   // smallest possible step size for integration;

  h[1]=1.0;
  for (j=1;j<=JMAX+1;j++) {
    s[j]=trapzd(x,y,hmin,xlength,j);
    if (j >= K) {
      lastss = ss;
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss) || fabs(lastss-ss)<1e-15) {
	//printf("number of attempts = %d\n",j);
	//printf("JMAX = %d\n",JMAX);
	*failure = false;
	//	printf("ss=%e\n",ss);
	return ss;
      }
    }
    h[j+1]=0.25*h[j];
    /*This is a key step: The factor is 0.25 even though the stepsize is decreased by only
0.5. This makes the extrapolation a polynomial in h^2 as allowed by equation (4.2.1),
not just a polynomial in h.*/
  }
  //  printf("Too many steps in routine qromb, retrying with %d steps.\n",(xlength-1)*2+1);
  //  printf("fabs(ss)*EPS = %e, fabs(dss) = %e, hmin = %e\n",EPS*fabs(ss),fabs(dss),hmin);
  //  printf("ss=%e\n",ss);
  *failure = true;
  return 0.0; // this return value doesn't matter, as failure signals that convergence failed.
}

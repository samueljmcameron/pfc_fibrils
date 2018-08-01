#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"

extern int mpt; //Defined in EvsR
extern double h, r[],K33,K24;
void difeq(int k, int k1, int k2, int jsf, int isl, int isf, int indexv[],
	   int ne, double **s, double **y)
/*Returns matrix s for solvde. See Numerical Recipes in C for details on what s is. */
{
  double temp1,temp2, shalf,chalf,sine,cosine;
  if (k == k1) { // BC at first point.
    s[2][2+indexv[1]] = 1.0;
    s[2][2+indexv[2]] = 0.0;
    s[2][jsf] = y[1][1];
  }
  else if (k > k2) { // BC at last point
    s[1][2+indexv[1]] = (-1.0*cos(2*y[1][mpt])+K24*2*cos(2*y[1][mpt]))/r[mpt];
    s[1][2+indexv[2]] = -1.0;
    s[1][jsf] = (1-y[2][mpt]-sin(2*y[1][mpt])/(2.0*r[mpt]))+K24*sin(2*y[1][mpt])/r[mpt];
  }
  else { // interior point
    temp1 = r[k]+r[k-1];
    temp2 = y[1][k]+y[1][k-1];
    sine = sin(temp2);
    cosine = cos(temp2);
    shalf = sin(temp2/2.0);
    chalf = cos(temp2/2.0);
    s[1][indexv[1]] = -1.0;
    s[1][indexv[2]] = -0.5*h;
    s[1][2+indexv[1]] = 1.0;
    s[1][2+indexv[2]] = -0.5*h;
    s[2][indexv[1]] = -h*(4*K33*shalf*sine*chalf/((temp1*temp1))\
			  +4*K33*(shalf*shalf)*cosine/((temp1*temp1))\
			  +2*sine*(1-sine/temp1)/temp1+2*(cosine*cosine)/(temp1*temp1));
    s[2][indexv[2]] = h/temp1-1;
    s[2][2+indexv[1]] = -h*(4*K33*shalf*sine*chalf/((temp1*temp1))\
			    +4*K33*(shalf*shalf)*cosine/((temp1*temp1)) \
			    +2*sine*(1-sine/temp1)/temp1+2*(cosine*cosine)/(temp1*temp1));
    s[2][2+indexv[2]] = h/temp1+1;
    s[1][jsf] = y[1][k]-y[1][k-1]-0.5*h*(y[2][k]+y[2][k-1]);
    s[2][jsf] = -h*((2*(1-(1.0/2.0)*y[2][k]-(1.0/2.0)*y[2][k-1]))/temp1\
		    +4*K33*(shalf*shalf)*sine/((temp1*temp1))\
		    -2*cosine*(1-sine/temp1)/temp1)+y[2][k]-y[2][k-1];
  }
}

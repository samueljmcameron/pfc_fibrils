#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"


void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,
	   struct params *p, double h, int mpt)
/*Returns matrix s for solvde. See Numerical Recipes in C for details on what s is. */
{
  double t1,t2, shalf,chalf,sine,cosine;
  double tmp;
  if (k == k1) { // BC at first point.
    s[2][3] = 1.0;
    s[2][4] = 0.0;
    s[2][jsf] = y[1][1];
  }
  else if (k > k2) { // BC at last point
    s[1][3] = -p->k24*cos(2*y[1][mpt])/r[mpt];
    s[1][4] = 1.0;
    s[1][jsf] = y[2][mpt]-1-p->k24*sin(2*y[1][mpt])/(2.0*r[mpt]);
  }
  else { // interior point
    tmp = 4*M_PI*M_PI/(p->d0*p->d0);
    t1 = r[k]+r[k-1];
    t2 = y[1][k]+y[1][k-1];
    sine = sin(t2);
    cosine = cos(t2);
    shalf = sin(t2/2.0);
    chalf = cos(t2/2.0);
    s[1][1] = -1.0;
    s[1][2] = -0.5*h;
    s[1][3] = 1.0;
    s[1][4] = -0.5*h;
    s[2][1] = -2*h/t1*(2*p->K33/t1*(sine*sine/2.0+shalf*shalf*cosine)
		       +sine*(1-sine/t1)+cosine*cosine/t1
		       +0.5*tmp*p->Lambda*p->delta*p->delta*t1
		       //		       *(1+sin(2*eta*L)/(2*eta*L))
		       *(tmp*shalf*shalf
			 /(chalf*chalf*chalf*chalf*chalf*chalf)
			 +0.5*(tmp/(chalf*chalf)-p->eta*p->eta)
			 *(1.0/(chalf*chalf)+3*shalf*shalf
			   /(chalf*chalf*chalf*chalf))));
    s[2][2] = h/t1-1;
    s[2][3] = -2*h/t1*(2*p->K33/t1*(sine*sine/2.0+shalf*shalf*cosine)
		       +sine*(1-sine/t1)+cosine*cosine/t1
		       +0.5*tmp*p->Lambda*p->delta*p->delta*t1
		       //		       *(1+sin(2*eta*L)/(2*eta*L))
		       *(tmp*shalf*shalf
			 /(chalf*chalf*chalf*chalf*chalf*chalf)
			 +0.5*(tmp/(chalf*chalf)-p->eta*p->eta)
			 *(1.0/(chalf*chalf)+3*shalf*shalf
			   /(chalf*chalf*chalf*chalf))));
    s[2][4] = h/t1+1;
    s[1][jsf] = y[1][k]-y[1][k-1]-0.5*h*(y[2][k]+y[2][k-1]);
    s[2][jsf] = (-2.0*h/t1*(1-0.5*(y[2][k]+y[2][k-1])
			    +2*p->K33*shalf*shalf*sine/t1
			    -cosine*(1-sine/t1)
			    +0.5*tmp*p->Lambda*p->delta*p->delta
			    //			    *(1+sin(2*eta*L)/(2*eta*L))
			    *(tmp/(chalf*chalf)-p->eta*p->eta)
			    *shalf*t1/(chalf*chalf*chalf))
		 +y[2][k]-y[2][k-1]);
  }
}

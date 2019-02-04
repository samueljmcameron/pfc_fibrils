#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"


void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   struct params *p,const double *x, double h)
/*Returns matrix s for solvde. See Numerical Recipes in C for details on what s is. */
{
  double t1,t2, shalf,chalf,sine,cosine;
  double tmp;
  if (k == k1) { // BC at first point.
    p->s[2][3] = 1.0;
    p->s[2][4] = 0.0;
    p->s[2][jsf] = p->y[1][1];
  }
  else if (k > k2) { // BC at last point
    p->s[1][3] = -p->k24*cos(2*p->y[1][p->mpt])/p->r[p->mpt];
    p->s[1][4] = 1.0;
    p->s[1][jsf] = p->y[2][p->mpt]-1-p->k24*sin(2*p->y[1][p->mpt])/(2.0*p->r[p->mpt]);
  }
  else { // interior point
    tmp = 4*M_PI*M_PI/(p->d0*p->d0);
    t1 = p->r[k]+p->r[k-1];
    t2 = p->y[1][k]+p->y[1][k-1];
    sine = sin(t2);
    cosine = cos(t2);
    shalf = sin(t2/2.0);
    chalf = cos(t2/2.0);
    p->s[1][1] = -1.0;
    p->s[1][2] = -0.5*h;
    p->s[1][3] = 1.0;
    p->s[1][4] = -0.5*h;
    p->s[2][1] = -2*h/t1*(2*p->K33/t1*(sine*sine/2.0+shalf*shalf*cosine)
			  +sine*(1-sine/t1)+cosine*cosine/t1
			  +0.5*tmp*p->Lambda*x[3]*x[3]*t1
			  *(tmp*shalf*shalf
			    /(chalf*chalf*chalf*chalf*chalf*chalf)
			    +0.5*(tmp/(chalf*chalf)-x[2]*x[2])
			    *(1.0/(chalf*chalf)+3*shalf*shalf
			      /(chalf*chalf*chalf*chalf))));
    p->s[2][2] = h/t1-1;
    p->s[2][3] = -2*h/t1*(2*p->K33/t1*(sine*sine/2.0+shalf*shalf*cosine)
			  +sine*(1-sine/t1)+cosine*cosine/t1
			  +0.5*tmp*p->Lambda*x[3]*x[3]*t1
			  *(tmp*shalf*shalf
			    /(chalf*chalf*chalf*chalf*chalf*chalf)
			    +0.5*(tmp/(chalf*chalf)-x[2]*x[2])
			    *(1.0/(chalf*chalf)+3*shalf*shalf
			      /(chalf*chalf*chalf*chalf))));
    p->s[2][4] = h/t1+1;
    p->s[1][jsf] = p->y[1][k]-p->y[1][k-1]-0.5*h*(p->y[2][k]+p->y[2][k-1]);
    p->s[2][jsf] = (-2.0*h/t1*(1-0.5*(p->y[2][k]+p->y[2][k-1])
			       +2*p->K33*shalf*shalf*sine/t1
			       -cosine*(1-sine/t1)
			       +0.5*tmp*p->Lambda*x[3]*x[3]
			       *(tmp/(chalf*chalf)-x[2]*x[2])
			       *shalf*t1/(chalf*chalf*chalf))
		    +p->y[2][k]-p->y[2][k-1]);
  }
}

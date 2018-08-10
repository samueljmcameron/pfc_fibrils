#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"


void solvde(int itmax, double conv, double slowc, double scalv[],
	    int ne, int nb, int m, double **y, double *r,
	    double ***c, double **s, double K33, double k24,
	    double Lambda,double d0,double L,double eta,double delta,
	    double h)
/*Driver routine for solution of two point boundary value problems by relaxation. itmax is the maximum number of iterations. conv is the convergence criterion (see text). slowc controls the fraction of corrections actually used after each iteration. scalv[1..ne] contains typical sizes for each dependent variable, used to weight errors. indexv[1..ne] lists the column ordering of variables used to construct the matrix s[1..ne][1..2*ne+1] of derivatives. (The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv.) The problem involves ne equations for ne adjustable dependent variables at each point. At the first mesh point there are nb boundary conditions. There are a total of m mesh points. y[1..ne][1..m] is the two-dimensional array that contains the initial guess for all the dependent variables at each mesh point. On each iteration, it is updated by the calculated correction. The arrays c[1..ne][1..ne-nb+1][1..m+1] and s supply dummy storage used by the relaxation code.*/
{

  void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
  void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	     int ne, double **s, double **y, double *r,double K33,
	     double k24, double Lambda, double d0,double L,
	     double eta, double delta, double h, int mpt);
  void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	     double ***c, double **s);
  void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	   int ic1, int jc1, int jcf, int kc, double ***c, double **s);
  int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
  int jc1,jcf,k,k1,k2,km,kp,nvars,*kmax;
  int i;
  double err,errj,fac,vmax,vz,*ermax;
  kmax=ivector(1,ne);
  ermax=vector(1,ne);
  k1=1;          //Set up row and column markers.
  k2=m;
  nvars=ne*m;
  j1=1;
  j2=nb;
  j3=nb+1;
  j4=ne;
  j5=j4+j1;
  j6=j4+j2;
  j7=j4+j3;
  j8=j4+j4;
  j9=j8+j1;
  ic1=1;
  ic2=ne-nb;
  ic3=ic2+1;
  ic4=ne;
  jc1=1;
  jcf=ic3;
  int index = 1;
  for (it=1;it<=itmax;it++) { //Primary iteration loop.
    k=k1; //Boundary conditions at first point.
    difeq(k,k1,k2,j9,ic3,ic4,ne,s,y,r,K33,k24,Lambda,d0,L,eta,delta,h,m);
    pinvs(ic3,ic4,j5,j9,jc1,k1,c,s);
    for (k=k1+1;k<=k2;k++) { //Finite difference equations at all point pairs.
      kp=k-1;
      difeq(k,k1,k2,j9,ic1,ic4,ne,s,y,r,K33,k24,Lambda,d0,L,eta,delta,h,m);
      red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
      pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
    }
    k=k2+1;// Final boundary conditions.
    difeq(k,k1,k2,j9,ic1,ic2,ne,s,y,r,K33,k24,Lambda,d0,L,eta,delta,h,m);
    red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s);
    pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s);
    bksub(ne,nb,jcf,k1,k2,c); //Backsubstitution.
    err=0.0;
    for (j=1;j<=ne;j++) { //Convergence check, accumulate average error
      errj=vmax=0.0;
      km=0;
      for (k=k1;k<=k2;k++) {// Find point with largest error, for each dependent variable
	vz=fabs(c[j][1][k]);
	if (vz > vmax) {
	  vmax=vz;
	  km=k;
	}
	errj += vz;
      }
      err += errj/scalv[j]; //Note weighting for each dependent variable.
      ermax[j]=c[j][1][km]/scalv[j];
      kmax[j]=km;
    }
    err /= nvars;
    fac=(err > slowc ? slowc/err : 1.0);
    //Reduce correction applied when error is large.
    for (j=1;j<=ne;j++) { //Apply corrections.
      for (k=k1;k<=k2;k++)
	y[j][k] -= fac*c[j][1][k];
    }
    //    printf("\n%8s %9s %9s\n","Iter.","Error","FAC"); //Summary of corrections
    //printf("%6d %12.12f %11.6f\n",it,err,fac);        //for this step.
    if (y[1][2]<0) {
      //      printf("bad!\n");
      //printf("psi(%lf) = %e\n",r[m],y[1][m]);
    }
    if (err < conv) { // Point with largest error for each variable can
                      // be monitored by writing out kmax and
                      // ermax.
      free_vector(ermax,1,ne);
      free_ivector(kmax,1,ne);
      return;
    }
  }
  nrerror("Too many iterations in solvde"); //Convergence failed.
}


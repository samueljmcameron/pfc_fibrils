#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"

double E2233b(double R,double *r,double *rf_,int mpt);
void compute_rf2233b(double K33, double k24, double Lambda,
		     double eta, double d0, double L, 
		     double *r, double **y, double *rf_,
		     int mpt);

void solvde(int itmax, double conv, double slowc, double scalv[],
	    int ne, int nb, int m, double **y, double *r,
	    double ***c, double **s, double K33, double k24,
	    double Lambda, double eta, double d0, double L,
	    double h,int mpt)
/*Driver routine for solution of two point boundary value problems by relaxation. itmax is the maximum number of iterations. conv is the convergence criterion (see text). slowc controls the fraction of corrections actually used after each iteration. scalv[1..ne] contains typical sizes for each dependent variable, used to weight errors. indexv[1..ne] lists the column ordering of variables used to construct the matrix s[1..ne][1..2*ne+1] of derivatives. (The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv.) The problem involves ne equations for ne adjustable dependent variables at each point. At the first mesh point there are nb boundary conditions. There are a total of m mesh points. y[1..ne][1..m] is the two-dimensional array that contains the initial guess for all the dependent variables at each mesh point. On each iteration, it is updated by the calculated correction. The arrays c[1..ne][1..ne-nb+1][1..m+1] and s supply dummy storage used by the relaxation code.*/
{

  void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
  void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	     int ne, double **s, double **y, double *r,double K33,
	     double k24, double Lambda, double eta, double d0,
	     double L,double h, int mpt);
  void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	     double ***c, double **s);
  void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	   int ic1, int jc1, int jcf, int kc, double ***c, double **s);
  int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
  int jc1,jcf,k,k1,k2,km,kp,nvars,*kmax;
  int i;
  double err,errj,fac,vmax,vz,*ermax;
  char f1[200],f2[200],f3[200],path[200];
  FILE *E_vs_it,*psi_vs_r,*rf_vs_r;
  double rf_[mpt+1];
  double Energy;

  snprintf(path,sizeof(path),"../Results/");
  snprintf(f1,sizeof(f1),"%sE_vs_it_k24_%1.4e_Lambda_%1.4e_R_%1.4e.txt",
	   path,k24,Lambda,r[mpt]);
  snprintf(f2,sizeof(f2),"%spsi_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e.txt",
	   path,k24,Lambda,r[mpt]);
  snprintf(f3,sizeof(f3),"%srf_vs_r_k24_%1.4e_Lambda_%1.4e_R_%1.4e.txt",
	   path,k24,Lambda,r[mpt]);

  E_vs_it = fopen(f1,"w");
  psi_vs_r = fopen(f2,"w");
  rf_vs_r = fopen(f3,"w");

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
  for (it=1;it<=itmax;it++) { //Primary iteration loop.
    k=k1; //Boundary conditions at first point.
    difeq(k,k1,k2,j9,ic3,ic4,ne,s,y,r,K33,k24,Lambda,eta,d0,L,h,mpt);
    pinvs(ic3,ic4,j5,j9,jc1,k1,c,s);
    for (k=k1+1;k<=k2;k++) { //Finite difference equations at all point pairs.
      kp=k-1;
      difeq(k,k1,k2,j9,ic1,ic4,ne,s,y,r,K33,k24,Lambda,eta,d0,L,h,mpt);
      red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
      pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
    }
    k=k2+1;// Final boundary conditions.
    difeq(k,k1,k2,j9,ic1,ic2,ne,s,y,r,K33,k24,Lambda,eta,d0,L,h,mpt);
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
    compute_rf2233b(K33,k24,Lambda,eta,d0,L,r,y,rf_,mpt);
    
    Energy = E2233b(r[mpt],r,rf_,mpt);
    printf("\nit = %d, E = %.12e\n",it,Energy);
    fprintf(E_vs_it,"%d\t%.12e\n",it,Energy);

    printf("\n%8s %9s %9s\n","Iter.","Error","FAC"); //Summary of corrections
    printf("%6d %12.12f %11.6f\n",it,err,fac);        //for this step.
    if (y[1][2]<0) {
      printf("bad!\n");
      printf("psi(R) = %e\n",y[1][mpt]);
    }
    if (err < conv) { // Point with largest error for each variable can
                      // be monitored by writing out kmax and
                      // ermax.
      for (i = 1; i <= mpt; i++) {
	fprintf(psi_vs_r,"%e\t%e\t%e\n",r[i],y[1][i],y[2][i]);
	fprintf(rf_vs_r,"%e\t%e\n",r[i],rf_[i]);
      }
      free_vector(ermax,1,ne);
      free_ivector(kmax,1,ne);
      return;
    }
  }
  fclose(E_vs_it);
  fclose(psi_vs_r);
  fclose(rf_vs_r);
  nrerror("Too many iterations in solvde"); //Convergence failed.
}

void compute_rf2233b(double K33, double k24, double Lambda,
		     double eta, double d0, double L, 
		     double *r, double **y, double *rf_,
		     int mpt)
{
  int i;
  double siny, sin2y;
  
  rf_[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r

    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    rf_[i] = r[i]*(-(y[2][i]+0.5*sin2y/r[i])
		   +0.5*(y[2][i]+0.5*sin2y/r[i])
		   *(y[2][i]+0.5*sin2y/r[i])
		   +0.5*K33*siny*siny*siny*siny/(r[i]*r[i])
		   +Lambda/4.0*(1+sin(2*eta*L)/(2*eta*L))
		   *(4*M_PI*M_PI/(d0*d0*cos(y[1][i])*cos(y[1][i]))
		     -eta*eta)
		   *(4*M_PI*M_PI/(d0*d0*cos(y[1][i])*cos(y[1][i]))
		     -eta*eta));
  }
  return;
}

double E2233b(double R,double *r,double *rf_,int mpt)
{
  return 2/(R*R)*qromb(r,rf_,mpt);
}

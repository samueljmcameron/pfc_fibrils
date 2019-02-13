#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define EPS 1.0e-14

void solvde_wrapper(double scalv[],struct params *p,double h,
		    bool ignore_first_y)
/*==============================================================================

  Purpose: Runs solvde up to three times. The first time using the y form which
  comes from previous calculations of y. If that fails, a linear guess for y 
  with a final y[1][mpt] = M_PI/4.0 If that fails, a psi(r) = a*sqrt(r) guess
  is used for y, with a = M_PI/(2.01*sqrt(R)). If that fails, give up and exit
  to system, while saving the initial guess for psi(r), as well as the current
  form of psi(r) which is not a solution to the ODE either (as convergence from
  the relaxation method in solvde failed).
  ------------------------------------------------------------------------------

  Parameters:

  ------------------------------------------------------------------------------

  Returns: Does not return

  ============================================================================*/
{

  bool solvde(double scalv[],struct params *p,double h,bool flag);

  void sqrtGuess(double *r, double **y, double initialSlope,double h,int mpt);

  void linearGuess(double *r, double **y, double initialSlope,double h,int mpt);

  void write_SOLVDEfailure(struct params *p);

  bool success_brent(double *psip01,double psip02,struct params *p,
		     double h);
  
  void constantGuess(double *r,double **y,double val,double h,int mpt);


  double slopeguess;
  double psip01 = 0;
  double psip02;
  bool flag = true;


  if (!solvde(scalv,p,h,flag) || ignore_first_y) {

    flag = true;//false;
    printf("solvde convergence failed when (R,eta,delta) = (%e,%e,%e).\n",p->R,p->eta,p->delta);
    //printf("Retrying using a hybrid shooting/relaxing method.\n");
    printf("Retrying using a constant y value\n");
    
    //  if (p->eta > 6.32) psip02 = M_PI/(0.01*p->R);
    //else psip02 = M_PI/(2.0*p->R);
    //success_brent(&psip01,psip02,p,h);
    constantGuess(p->r,p->y,0.3,h,p->mpt);

  } else return;
  
  if (!solvde(scalv,p,h,flag)) {

    flag = true;
    printf("Retrying with a linear guess and a final twist angle value of "
	   "pi/4.\n");

    slopeguess = M_PI/(4.0*p->R);
    
    linearGuess(p->r,p->y,slopeguess,h,p->mpt);

  } else return;
  if (!solvde(scalv,p,h,flag)) {

    flag = true;
    printf("Retrying with a sqrt(r) guess and a final twist angle value of "
	   "pi/2.01.\n");

    slopeguess = M_PI/(2.01*sqrt(p->R));

    sqrtGuess(p->r,p->y,slopeguess,h,p->mpt);

  } else return;


  if (!solvde(scalv,p,h,flag)) {
    
    // save form of y when solvde failed, rf_fib, and exit.


    write_SOLVDEfailure(p);
  } else return;

  
}

bool solvde(double scalv[],struct params *p,double h,bool flag)
/*==============================================================================
  Driver routine for solution of two point boundary value problems
  by relaxation. itmax is the maximum number of iterations. conv
  is the convergence criterion (see text). slowc controls the fraction
  of corrections actually used after each iteration. scalv[1..ne]
  contains typical sizes for each dependent variable, used to weight
  errors. indexv[1..ne] lists the column ordering of variables used to
  construct the matrix s[1..ne][1..2*ne+1] of derivatives. (The nb
  boundary conditions at the first mesh point must contain some 
  dependence on the first nb variables listed in indexv.) The problem
  involves ne equations for ne adjustable dependent variables at each
  point. At the first mesh point there are nb boundary conditions. 
  There are a total of m mesh points. y[1..ne][1..m] is the two-
  dimensional array that contains the initial guess for all the 
  dependent variables at each mesh point. On each iteration, it is 
  updated by the calculated correction. The arrays
  c[1..ne][1..ne-nb+1][1..m+1] and s supply dummy storage used by the 
  relaxation code.*/
{

  void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
  void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	     struct params *p,double h);
  bool pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	     double **s);
  void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	   int ic1, int jc1, int jcf, int kc, double ***c, double **s);


  int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
  int jc1,jcf,k,k1,k2,km,kp,nvars,*kmax;
  int i;
  double err,errj,fac,vmax,vz,*ermax;
  double slowc = 1.0;

  kmax=ivector(1,NE);
  ermax=vector(1,NE);
  k1=1;          //Set up row and column markers.
  k2=p->mpt;
  nvars=NE*p->mpt;
  j1=1;
  j2=NB;
  j3=NB+1;
  j4=NE;
  j5=j4+j1;
  j6=j4+j2;
  j7=j4+j3;
  j8=j4+j4;
  j9=j8+j1;
  ic1=1;
  ic2=NE-NB;
  ic3=ic2+1;
  ic4=NE;
  jc1=1;
  jcf=ic3;

  for (it=1;it<=ITMAX;it++) { //Primary iteration loop.
    k=k1; //Boundary conditions at first point.
    difeq(k,k1,k2,j9,ic3,ic4,p,h);
    //    if (isnan(y[1][k])) printf("NAN at first BC!\n");
    if (!pinvs(ic3,ic4,j5,j9,jc1,k1,p->c,p->s)) {
      printf("R = %e\n",p->R);
      printf("failed at first BC!\n");
      return false;
    }
    for (k=k1+1;k<=k2;k++) { //Finite difference equations at all point pairs.
      kp=k-1;
      difeq(k,k1,k2,j9,ic1,ic4,p,h);
      //      if (isnan(y[1][k])) printf("NAN at k = %d!\n",k);
      red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,p->c,p->s);
      if (!pinvs(ic1,ic4,j3,j9,jc1,k,p->c,p->s)) {
	printf("R = %e\n",p->R);
	printf("failed at point k = %d in finite differences\n",k);
	return false;
      }
    }
    k=k2+1;// Final boundary conditions.
    difeq(k,k1,k2,j9,ic1,ic2,p,h);
    //    if (isnan(y[1][k])) printf("NAN at last BC!\n");
    red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,p->c,p->s);
    if (!pinvs(ic1,ic2,j7,j9,jcf,k2+1,p->c,p->s)) {
      printf("R = %e\n",p->R);
      printf("failed at last BC!\n");
      return false;
    }
    bksub(NE,NB,jcf,k1,k2,p->c); //Backsubstitution.
    err=0.0;
    for (j=1;j<=NE;j++) { //Convergence check, accumulate average error
      errj=vmax=0.0;
      km=0;
      for (k=k1;k<=k2;k++) {// Find point with largest error, for each dependent variable
	vz=fabs(p->c[j][1][k]);
	if (vz > vmax) {
	  vmax=vz;
	  km=k;
	}
	errj += vz;
      }
      err += errj/scalv[j]; //Note weighting for each dependent variable.
      ermax[j]=p->c[j][1][km]/scalv[j];
      kmax[j]=km;
    }
    err /= nvars;
    fac=(err > slowc ? slowc/err : 1.0);
    //Reduce correction applied when error is large.
    for (j=1;j<=NE;j++) { //Apply corrections.
      for (k=k1;k<=k2;k++) {
	p->y[j][k] -= fac*p->c[j][1][k];
	//if (isnan(p->y[j][k])) printf("NAN!\n");
      }
    }
    if (err < CONV_ODE) { // Point with largest error for each variable can
                      // be monitored by writing out kmax and
                      // ermax.

      if (flag) {
	if (p->y[1][2]<=-1e-15) {
	  printf("likely a maximizing (vs minimizing) solution for psi(r), "
		 "as y[1][2]=%e  which usually means high energy.\n",p->y[1][2]);
	  return false;
	}
      }

      free_vector(ermax,1,NE);
      free_ivector(kmax,1,NE);
      return true;
    }
  }
  printf("Too many iterations in solvde, err = %e.\n",err); //Convergence failed.
  return false;
}

void constantGuess(double *r,double **y,double val,double h,int mpt)
{
  int k;
  for (k = 1; k <= mpt; k++) {
    r[k] = (k-1)*h;
    y[1][k] = val;
    y[2][k] = 0.0;
  }
  return;
}

void sqrtGuess(double *r, double **y, double initialSlope,double h,int mpt)
/*==============================================================================

  Purpose:

  ------------------------------------------------------------------------------

  Parameters:

  ------------------------------------------------------------------------------

  Returns: Does not return anything explicitly, just modifies the y array.

  ============================================================================*/

{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess!
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*sqrt(r[k]); // y1 is psi
    y[2][k] = initialSlope/(2.0*sqrt(r[k]+0.01)); // y2 is psi'
  }
  return;
}

void linearGuess(double *r, double **y, double initialSlope,double h,int mpt)
/*==============================================================================

  Purpose:

  ------------------------------------------------------------------------------

  Parameters:

  ------------------------------------------------------------------------------

  Returns: Does not return anything explicitly, just modifies the y array.

  ============================================================================*/

{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*r[k]; // y1 is psi
    y[2][k] = initialSlope; // y2 is psi'
  }
  return;
} 
  
void write_SOLVDEfailure(struct params *p)
/*==============================================================================
  
  Purpose: This function saves the current forms of r and psi(r) (which has
  not relaxed to the correct psi(r) form), as well as the form used as the 
  initial guess for psi(r) which is inputted to the solvde_wrapper. This
  function is only called if solvde_wrapper fails to calculate psi(r). 

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ------------------------------------------------------------------------------

  Returns: Exits with exit status exit(1);

  ============================================================================*/

{

  void make_f_err(char *f_err,char *err_type,int f_err_size,struct params *p);

  int i;
  FILE *broken1,*broken2;
  int f_err_size = 200;
  char f_err1[f_err_size];
  char f_err2[f_err_size];

  printf("failed to solve ODE when (R,eta,delta) = (%e,%e,%e).\n",p->R,p->eta,p->delta);
  printf("saving current psi(r) shape (from failed solvde call), as well as the"
	 " shape of the initial guess of psi(r) (from previous call of E_calc),"
	 " and exiting to system.\n");

  make_f_err(f_err1,"SOLVDE_FAIL",f_err_size,p);
  broken1 = fopen(f_err1,"w");

  for (i = 1; i<=p->mpt; i++) {
    fprintf(broken1,"%.8e\t%.8e\t%.8e\n",p->r[i],p->y[1][i],p->y[2][i]);
  }
  fclose(broken1);

  make_f_err(f_err2,"SOLVDE_INITGUESS",f_err_size,p);
  broken2 = fopen(f_err2,"w");

  for (i = 1; i<=p->mpt; i++) {
    fprintf(broken2,"%.8e\t%.8e\t%.8e\n",p->r[i],p->y_cp[1][i],p->y_cp[2][i]);
  }
  fclose(broken2);

  exit(1);
  return;
}


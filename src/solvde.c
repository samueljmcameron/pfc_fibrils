#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"



void solvde_wrapper(int itmax, double conv, double scalv[],struct arr_ns *ns,
		    int mpt,double *r,double **y,double **y_guess,double ***c,
		    double **s,struct params *p,double *x,double h)
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

  bool solvde(int itmax, double conv, double scalv[],struct arr_ns *ns, int m,
	      double *r, double **y,double ***c, double **s,struct params *p,
	      double *x,double h);

  void sqrtGuess(double *r, double **y, double initialSlope,double h,int mpt);

  void write_SOLVDEfailure(double *r,double **y,double **y_guess,int mpt,
			   struct params p,double *x);


  double slopeguess;



  if (!solvde(itmax,conv,scalv,ns,mpt,r,y,c,s,p,x,h)) {
    printf("solvde convergence failed, trying one more time with a "
	   "linear guess and a final twist angle value of pi/4.\n");

    slopeguess = M_PI/(4.0*x[1]);
    
    linearGuess(r,y,slopeguess,h,mpt);

  } else return;
  if (!solvde(itmax,conv,scalv,ns,mpt,r,y,c,s,p,x,h)) {

    
    printf("solvde convergence failed, trying one more time with a "
	   "sqrt(r) guess and a final twist angle value of pi/(2.01).\n");

    slopeguess = M_PI/(2.01*sqrt(x[1]));

    sqrtGuess(r,y,slopeguess,h,mpt);

  } else return;
  if (!solvde(itmax,conv,scalv,ns,mpt,r,y,c,s,p,x,h)) {
    
    // save form of y when solvde failed, rf_fib, and exit.


    write_SOLVDEfailure(r,y,y_guess,mpt,*p,x);
    


  } else return;
}

bool solvde(int itmax, double conv, double scalv[],struct arr_ns *ns, int m,
	    double *r, double **y,double ***c, double **s,struct params *p,
	    double *x,double h)
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
  /*
  void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
  bool pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	     double ***c, double **s);
  void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	   int ic1, int jc1, int jcf, int kc, double ***c, double **s);
  */
  int ne = ns->ne;
  int nb = ns->nb;
  int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
  int jc1,jcf,k,k1,k2,km,kp,nvars,*kmax;
  int i;
  double err,errj,fac,vmax,vz,*ermax;
  double slowc = 1.0;

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
    difeq(k,k1,k2,j9,ic3,ic4,ne,s,y,r,p,x,h,m);
    //    if (isnan(y[1][k])) printf("NAN at first BC!\n");
    if (!pinvs(ic3,ic4,j5,j9,jc1,k1,c,s)) {
      printf("R = %e\n",x[1]);
      printf("failed at first BC!\n");
      return false;
    }
    for (k=k1+1;k<=k2;k++) { //Finite difference equations at all point pairs.
      kp=k-1;
      difeq(k,k1,k2,j9,ic1,ic4,ne,s,y,r,p,x,h,m);
      //      if (isnan(y[1][k])) printf("NAN at k = %d!\n",k);
      red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
      if (!pinvs(ic1,ic4,j3,j9,jc1,k,c,s)) {
	printf("R = %e\n",x[1]);
	printf("failed at point k = %d in finite differences\n",k);
	return false;
      }
    }
    k=k2+1;// Final boundary conditions.
    difeq(k,k1,k2,j9,ic1,ic2,ne,s,y,r,p,x,h,m);
    //    if (isnan(y[1][k])) printf("NAN at last BC!\n");
    red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s);
    if (!pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)) {
      printf("R = %e\n",x[1]);
      printf("failed at last BC!\n");
      return false;
    }
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
      for (k=k1;k<=k2;k++) {
	y[j][k] -= fac*c[j][1][k];
	if (isnan(y[j][k])) printf("NAN!\n");
      }
    }
    //    printf("\n%8s %9s %9s\n","Iter.","Error","FAC"); //Summary of corrections
    //printf("%6d %12.12f %11.6f\n",it,err,fac);        //for this step.
    if (err < conv) { // Point with largest error for each variable can
                      // be monitored by writing out kmax and
                      // ermax.

      if (y[1][2]<=-1e-15) {
	printf("likely a maximizing (vs minimizing) solution for psi(r), "
	       "as y[1][2]=%e  which usually means high energy.\n",y[1][2]);
	return false;
      }

      free_vector(ermax,1,ne);
      free_ivector(kmax,1,ne);
      return true;
    }
  }
  printf("Too many iterations in solvde, err = %e.\n",err); //Convergence failed.
  return false;
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
  
void write_SOLVDEfailure(double *r,double **y,double **y_guess,int mpt,
			 struct params p,double *x)
/*==============================================================================
  
  Purpose: This function saves the current forms of r and psi(r) (which has
  not relaxed to the correct psi(r) form), as well as the form used as the 
  initial guess for psi(r) which is inputted to the solvde_wrapper. This
  function is only called if solvde_wrapper fails to calculate psi(r). 

  ------------------------------------------------------------------------------

  Parameters:

  r[1..max_mpt] -- This vector holds the grid points for psi(r). Only the 
  first mpt values, with r[mpt] = R (= x[1]) are used, but the remaining grid
  values are there in case interpolation needs to occur.

  y[1..2][1..max_mpt] -- This 2d matrix holds psi(r) in y[1][:] and dpsi/dr in
  y[2][:]. Again, only the first mpt values are used until interpolation is 
  necessary.

  y_guess -- Array with the initial guess of psi(r) and dpsi/dr.

  mpt -- The number of grid points in the r and psi(r) discretization.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Exits with exit status exit(1);

  ============================================================================*/

{



  int i;
  FILE *broken1,*broken2;
  int f_err_size = 200;
  char f_err1[f_err_size];
  char f_err2[f_err_size];

  printf("failed to solve ODE when x = (%e,%e,%e).\n",x[1],x[2],x[3]);
  printf("saving current psi(r) shape (from failed solvde call), as well as the"
	 " shape of the initial guess of psi(r) (from previous call of E_calc),"
	 " and exiting to system.\n");

  make_f_err(f_err1,"SOLVDE_FAIL",f_err_size,p,x);
  broken1 = fopen(f_err1,"w");

  for (i = 1; i<=mpt; i++) {
    fprintf(broken1,"%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i]);
  }
  fclose(broken1);

  make_f_err(f_err2,"SOLVDE_INITGUESS",f_err_size,p,x);
  broken2 = fopen(f_err2,"w");

  for (i = 1; i<=mpt; i++) {
    fprintf(broken2,"%.8e\t%.8e\t%.8e\n",r[i],y_guess[1][i],y_guess[2][i]);
  }
  fclose(broken2);

  exit(1);
  return;
}


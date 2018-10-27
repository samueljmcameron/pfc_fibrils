/* gradient descent derivatives for energy function */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]

// free energy density which includes only the terms which require integration //
// (i.e. the q,K22,K33, and Lambda terms) //

double f2233b1_r(struct params *p,double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp);


// calculation of integrand array to be put into qromb for integration. //

void compute_rf2233b1(struct params *p,double *x,double *r,
		      double **y, double *rf_fib,int mpt);

// calculation of energy, given that the integral has been calculated. //

double pre_integral_E(struct params *p,double *x,double *r,double **y,
		      double integration_2233b1,int mpt);

// calculation of the energy, returns true if integral was successfully
// calculated.

bool integrate_E(double *E,struct params *p,double *x,double *r,
		 double **y,double *rf_fib,int mpt);


void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x);

void write_QROMBfailure(double *r, double **y,double *rf_fib,
			int mpt,struct params p,double *x);

void write_SOLVDEfailure(double *r, double **y,double *r_cp, double **y_cp,
			 int mpt,int last_mpt,struct params p,double *x);


double f2233b1_r(struct params *p,double *x,double ri,double sin_yi,
		 double sin_2yi,double cos_yi,double yi_p,
		 double tmp)
{
  double ans;
  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri)+p->Lambda*x[3]*x[3]/4.0
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2])
	 *(tmp/(cos_yi*cos_yi)-x[2]*x[2]));

  return ans;
}

void compute_rf2233b1(struct params *p,double *x,double *r,
		      double **y, double *rf_fib,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  rf_fib[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    rf_fib[i] = r[i]*f2233b1_r(p,x,r[i],siny,sin2y,cosy,y[2][i],tmp);
  }
  return;
}

double pre_integral_E(struct params *p,double *x,double *r,double **y,
		      double integration_2233b1,int mpt)
{

  double E;



  E = 2.0/(x[1]*x[1])*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = (E+x[3]*x[3]*p->omega*0.5
       *(0.75*x[3]*x[3]-1
	 //	 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
	 //	 +delta*delta*sin(4*eta*L)/(16*eta*L)
	 ));
  // add surface term tension terms
  E = E+0.5+1.0/x[1]*(-(1+p->k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/x[1]+2.0*p->gamma_s);  
  //  E = E+2*gamma_t/L;

  return E;
}



bool integrate_E(double *E,struct params *p,double *x,double *r,
		 doxuble **y,double *rf_fib,int mpt)
// computes E(R) given psi(r) (i.e. y) values //
  
{
  
  bool failure = true;
  double integration_2233b1;
  double tol0 = 1e-14;
  double tol2233b1;

  // tolerances ("tol...") are to determine how large the energy or 
  // derivative terms are without the integrals (setting them to 0).
  // The reason for this is if the integral calculation performed
  // by qromb has not converged (usually because the integral error
  // is so small that round-off error becomes an issue), if the
  // size of the error term is so small relative to the rest of the
  // function that it doesn't change the functions value (up to
  // some tolerance tolblabla), then the lack of convergence can just
  // be ignored.

  
  tol2233b1 = fabs(pre_integral_E(p,r,y,0,mpt)*x[1]*x[1]/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,x,r,y,rf_fib,mpt);
  integration_2233b1 = qromb(r,rf_fib,mpt,tol2233b1,&failure);
  
  if (failure) {
    printf("failed to integrate at x = (%e,%e,%e)\n"
	   x[1],x[2],x[3],mpt);
    return false;
  }

  *E = pre_integral_E(p,x,r,y,integration_2233b1,mpt);

  return true;
}

double E_calc(struct params *p,double *x,double ***c,double **s,
	      double **y,double *r,double *rf_fib,double **y_cp,
	      double *r_cp,double conv,int itmax,int *mpt,
	      int last_mpt,struct arr_ns *ns,int max_size)
{
  double h;
  double E;

  double scalv[2+1];

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  while ((*mpt) <= max_size) {

    h = x[1]/((*mpt)-1);    // compute stepsize in r[1..mpt] 

    if (need_to_interpolate(*mpt,last_mpt)) {

      printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	     x[1],x[2],x[3]);

      copy_2_arrays(r,y,r_cp,y_cp,last_mpt); // copy arrays r and y into r_cp and y_cp
      propagate_r(r,h,*mpt);
      interpolate_array(r,y,r_cp,y_cp,*mpt); // interpolate old y values (now stored in y_cp)


      last_mpt = *mpt;
    }

    else propagate_r(r,h,(*mpt));
    
    solvde_wrapper(itmax,conv,scalv,ns,*mpt,y,r,c,s,p,x,h);

    if(integrate_E(*E,p,x,r,y,rf_fib,mpt)) return E;

    
    (*mpt) = ((*mpt)-1)*2+1;
    
  }

  // if it makes it this far, we did not successfully compute E(R)

  write_QROMBfailure(r,y,rf_fib,*mpt,*p,x); // save psi(r), rf_fib(r), and exit

  // never get here.

  return 0;
}

void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x)
{
  snprintf(f_err,f_err_size,"data/%s_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",err_type,
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,x[1],x[2],x[3],
	   p.gamma_s);
  return;
}

void write_QROMBfailure(double *r, double **y,double *rf_fib,
			int mpt,struct params p,double *x)
{
  int i;
  FILE *broken;
  int f_err_size = 200;
  char f_err[f_err_size];

  make_f_err(f_err,"QROMB",f_err_size,p,x);

  printf("failed to integrate with qromb at x = (%e,%e,%e).\n",x[1],x[2],x[3]);
  printf("saving psi(r) shape, rf_fib(r), and exiting to system.\n");

  broken = fopen(f_err,"w");


  for (i = 1; i<=rlength; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i],rf_fib[i]);
  }
  fclose(broken);
  exit(1);
  return;
}

void write_SOLVDEfailure(double *r, double **y,double *r_cp, double **y_cp,
			 int mpt,int last_mpt,struct params p,double *x)
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
    fprintf(broken,"%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i]);
  }
  fclose(broken1);

  make_f_err(f_err2,"SOLVDE_INITGUESS",f_err_size,p,x);
  broken2 = fopen(f_err2,"w");

  for (i = 1; i<=last_mpt; i++) {
    fprintf(broken2,"%.8e\t%.8e\t%.8e\n",r_cp[i],y_cp[1][i],y_cp[2][i]);
  }
  fclose(broken2);

  exit(1);
  return;
}


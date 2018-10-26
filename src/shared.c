#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"
#include <time.h>

void sqrtGuess(double *r, double **y, double initialSlope,
	       double h,int mpt);

bool need_to_interpolate(int mpt, int last_mpt);

void single_calc(double *E,double *dEdx,double *x,struct params *p,
		 double *x,double ***c,double **s,double **y,double *r,
		 double *rf_fib,double **y_cp,double *r_cp,double *hessian,
		 double conv,int itmax,int *mpt,int last_mpt,
		 struct arr_ns *ns,int max_size)
{
  double h;

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

    if(calc_E(E,p,x,r,y,rf_fib,mpt)) return;

    
    (*mpt) = ((*mpt)-1)*2+1;
    
  }

  // if it makes it this far, we did not successfully compute E(R)

  write_QROMBfailure(r,y,rf_fib,*mpt,*p,x); // save psi(r), rf_fib(r), and exit

}

bool need_to_interpolate(int mpt, int last_mpt)
{
  return mpt != last_mpt;
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
	 " shape of the initial guess of psi(r) (from previous call of single_calc),"
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

void save_psi(FILE *psi,double *r, double **y,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(psi,"%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i]);
  }
  printf("psi(R) = %1.2e\n",y[1][mpt]);
  return;
}

void save_energydensity(FILE *energydensity,double *r, double *rf_fib,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(energydensity,"%10.8e\t%10.8e\n",r[i],rf_fib[i]);
  }
  return;
}

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable)
{
  fprintf(energy,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,E,
	  derivative,observable);
  return;
}

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int mpt)
{
  int i;
  double dy;

  for (i = 1; i <=mpt; i++) quick_interp(r_cp,y_cp,r[i],y,i);
  printf("done interpolating.\n");

  return;
}


void assign_ns(struct arr_ns *ns)
{
  int ne = 2;
  int nb = 1;
  ns->ne = ne;
  ns->nb = nb;
  ns->nsi = ne;
  ns->nsj = 2*ne+1;
  ns->nyj = ne;
  ns->nci = ne;
  ns->ncj = ne-nb+1;
  return;
}


void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_fib,int mpt,struct arr_ns ns)
{
  *y = matrix(1,ns.nyj,1,mpt);
  *s = matrix(1,ns.nsi,1,ns.nsj);
  *c = f3tensor(1,ns.nci,1,ns.ncj,1,mpt+1);
  *r = vector(1,mpt);
  *rf_fib = vector(1,mpt);
  return;
}

void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_fib,int mpt,struct arr_ns ns)
{

  free_f3tensor(*c,1,ns.nci,1,ns.ncj,1,mpt+1);
  free_matrix(*s,1,ns.nsi,1,ns.nsj);
  free_matrix(*y,1,ns.nyj,1,mpt);
  free_vector(*r,1,mpt);
  free_vector(*rf_fib,1,mpt);
  return;
}

void sqrtGuess(double *r, double **y, double initialSlope,
		 double h,int mpt)
{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess!
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*sqrt(r[k]); // y1 is psi
    y[2][k] = initialSlope/(2.0*sqrt(r[k]+0.01)); // y2 is psi'!!!!!!!!!!!!!!!!!!
  }
  return;
}


void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt)
{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess
    //    printf("r[%d] = %lf\n",k,r[k]);
    r[k] = (k-1)*h;
    //    printf("r[%d] = %lf\n\n",k,r[k]);
    y[1][k] = initialSlope*r[k]; // y1 is psi
    y[2][k] = initialSlope; // y2 is psi'!!!!!!!!!!!!!!!!!!
  }
  return;
}

void propagate_r(double *r, double h,int mpt)
{
  int k;
  for (k=1;k <=mpt; k++) r[k] = (k-1)*h; // only change r since psi, psi' are stored from last loop
  return;
}

int index(char scan_what[])
{
  if (strcmp(scan_what,"R")==) {
    printf("R!\n");
    return 1;
  } else if (strcmp(scan_what,"eta")==0) {
    printf("eta!\n");
    return 2;
  } else if (strcmp(scan_what,"delta")==0) {
    printf("delta!\n");
    return 3;
  } else {
    printf("Need either R, eta, or delta as argv[n] input."
	   "Exiting to system.\n");
    exit(1);
  }
  return 0; // never get here
}


void quick_interp(double *xcp,double **ycp,double x,double **y,int i)
// Given two points of ycp[1][:] (y11cp,x1cp) and y(12cp,x2cp), //
// and two points of ycp[2][:] (y21cp,x1cp) and (y22cp,x2cp),   //
// interpolate to determine the value of y between the two

{
  if (i%2 != 0) {
    y[1][i] = ycp[1][(i-1)/2+1];
    y[2][i] = ycp[2][(i-1)/2+1];
  }
  else {
    double slope_1, slope_2;
    slope_1 = (ycp[1][i/2+1]-ycp[1][i/2])/(xcp[i/2+1]-xcp[i/2]);
    slope_2 = (ycp[2][i/2+1]-ycp[2][i/2])/(xcp[i/2+1]-xcp[i/2]);

    y[1][i] = slope_1*(x-xcp[i/2])+ycp[1][i/2];
    y[2][i] = slope_2*(x-xcp[i/2])+ycp[2][i/2];
  }

  return;
}


int sign(double x) 
{
  return (x > 0) - (x < 0);
}

void arr_cp(double *acp, double *a,int length)
{
  int i;
  for (i = 1; i <= length; i++) acp[i] = a[i];
  return;
}

double vector_norm(double *a,int length)
{
  int i;
  double sum = 0;
  for (i = 1; i <= length; i++) sum += a[i]*a[i];

  return sum;
}

void array_constant(double constant,double *a,int length)
{
  int i;
  for (i = 1; i <= length; i++) a[i] = constant;
  return;
}

bool non_zero_array(double *dEdx,double conv)
{
  return (fabs(dEdx[1]) > conv || fabs(dEdx[2]) > conv
	  || fabs(dEdx[3]) > conv);
}

void copy_2_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_mpt)
{
  int i;

  for (i = 1; i <=last_mpt; i++) {
    y_cp[1][i] = y[1][i];
    y_cp[2][i] = y[2][i];
    r_cp[i] = r[i];
  }

  return;
}

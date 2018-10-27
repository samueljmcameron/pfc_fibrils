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


bool need_to_interpolate(int mpt, int last_mpt)
{
  return mpt != last_mpt;
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

void allocate_vectors(double x_size,double **x, double **dEdx, double **lastdEdx,
		      double **direction, double **hessian, double **x_cp,
		      double **E_p, double **E_m, double **E_pij,
		      double **E_mij)
{

  // these vectors are meaningful
  x = vector(1,x_size); // x = (R,eta,delta)
  dEdx = vector(1,x_size); // dEdx = grad(E)
  lastdEdx = vector(1,x_size); // last grad(E)
  direction = vector(1,x_size); // descent direction to be taken
  hessian = vector(1,x_size*x_size); // flattened hessian matrix

  // these vectors are just dummy vectors to use in calculating
  // derivatives
  x_cp = vector(1,x_size);
  E_p = vector(1,x_size);
  E_m = vector(1,x_size);
  E_pij = vector(1,x_size);
  E_mij = vector(1,x_size);


  return;
}

void free_vectors(double x_size,double **x, double **dEdx, double **lastdEdx,
		  double **direction, double **hessian, double **x_cp,
		  double **E_p, double **E_m, double **E_pij,
		  double **E_mij)
{
  free_vector(x,1,x_size);
  free_vector(dEdx,1,x_size);
  free_vector(lastdEdx,1,x_size);
  free_vector(direction,1,x_size);
  free_vector(hessian,1,x_size*x_size);

  free_vector(x_cp,1,x_size);
  free_vector(E_p,1,x_size);
  free_vector(E_m,1,x_size);
  free_vector(E_pij,1,x_size);
  free_vector(E_mij,1,x_size);
  return;
}

void allocate_matrices(struct arr_ns ns,double ****c,double ***s,
		       double ***y,double **r,double ***y_cp,
		       double **r_cp,double **rf_fib,int mpt)
{
  *y = matrix(1,ns.nyj,1,mpt);
  *s = matrix(1,ns.nsi,1,ns.nsj);
  *c = f3tensor(1,ns.nci,1,ns.ncj,1,mpt+1);
  *r = vector(1,mpt);
  *y_cp = matrix(1,ns.nyj,1,mpt);
  *r_cp = vector(1,mpt);
  *rf_fib = vector(1,mpt);
  return;
}


void free_matrices(struct arr_ns ns,double ****c,double ***s,
		   double ***y,double **r,double ***y_cp,
		   double **r_cp,double **rf_fib,int mpt)
{

  free_f3tensor(*c,1,ns.nci,1,ns.ncj,1,mpt+1);
  free_matrix(*s,1,ns.nsi,1,ns.nsj);
  free_matrix(*y,1,ns.nyj,1,mpt);
  free_vector(*r,1,mpt);
  free_matrix(*y_cp,1,ns.nyj,1,mpt);
  free_vector(*r_cp,1,mpt);
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

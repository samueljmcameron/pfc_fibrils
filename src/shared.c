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

void single_calc(double *E,double *dEdx,double *x,struct params *p,
		 double ****c,double ***s,double ***y,double **r,
		 double **rf_fib,double **y_cp,double *r_cp,double *hessian,
		 double conv,int itmax,int *npoints,int last_npoints,
		 struct arr_ns *ns,int max_size)
{
  double h;
  bool successful_qromb=false;
  bool successful_solvde;
  double slopeguess = M_PI/(4.0*x[1]);
  double slowc = 1.0;
  double scalv[2+1];

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  do {
    h = x[1]/((*npoints)-1);    

    if ((*npoints) != last_npoints) {

      // if the first calculation for the current variable values is unsuccessful
      printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	     x[1],x[2],x[3]);

      resize_and_interp(h,c,s,y,r,rf_fib,y_cp,r_cp,*npoints,
			last_npoints,ns);

      printf("done interpolating.\n");
      last_npoints = *npoints;
    }
    
    // if last loop was successful with last_xpoint sized arrays, then just
    // use last y values as initial guess
    else propagate_r(*r,h,(*npoints));
    

    successful_solvde = solvde(itmax,conv,slowc,scalv,
			       ns,(*npoints),*y,*r,*c,*s,p,h); // relax to compute psi,
    //                                                            psi' curves, 
      


    if (!successful_solvde) {
      printf("solvde convergence failed, trying one more time with a "
	     "linear guess and a final twist angle value of pi/4.\n");

      linearGuess(*r,*y,slopeguess,h,(*npoints));
      successful_solvde = solvde(itmax,conv,slowc,scalv,
				 ns,(*npoints),*y,*r,*c,*s,p,h); // relax to compute psi,
      //                                                            psi' curves, 
    }
    if (!successful_solvde) {
      slopeguess = M_PI/(2.01*sqrt(x[1]));
      printf("solvde convergence failed, trying one more time with a "
	     "sqrt(r) guess and a final twist angle value of pi/(2.01).\n");
      sqrtGuess(*r,*y,slopeguess,h,(*npoints));
      successful_solvde = solvde(itmax,conv,slowc,scalv,
				 ns,(*npoints),*y,*r,*c,*s,p,h); // relax to compute psi,
      //                                                            psi' curves, 
    }
    if (!successful_solvde) {
      write_failure("SOLVDE",*r,*y,*rf_fib,*npoints,*p);
      return;
    }


    if (calc_Hess) {
      // calculate energy, derivatives (see energy.c for code)
      successful_qromb = energy_prop_with_hessian(E,dEdx,p,*r,*y,*rf_fib,
						  (*npoints),hessian);
    }
    else {
      // calculate energy, derivatives (see energy.c for code)
      successful_qromb = energy_properties(E,dEdx,p,*r,*y,*rf_fib,
					   (*npoints));
    }

    
    
    (*npoints) = ((*npoints)-1)*2+1;
    
  }
  while (!successful_qromb && (*npoints) <= max_size);



  (*npoints) = ((*npoints)-1)/2+1;


  if (!successful_qromb) { // writes the psi(r), r*f, etc, and exits
    write_failure("QROMB",*r,*y,*rf_fib,*npoints,*p);
  }


  return;
}


void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p)
{
  snprintf(f_err,f_err_size,"data/%s_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",err_type,
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,
	   p.delta,p.gamma_s);
  return;
}

void write_failure(char *err_type,double *r, double **y,double *rf_fib,int rlength,struct params p)
{
  int i;
  FILE *broken;
  int f_err_size = 200;
  char f_err[f_err_size];

  make_f_err(f_err,err_type,f_err_size,p);

  if (strcmp(err_type,"QROMB")==0) {
      printf("failed to integrate at npoints = %d, with R = %e.\n",rlength,r[rlength]);
      printf("saving psi(r) shape, and exiting to system.\n");
    }
  broken = fopen(f_err,"w");


  for (i = 1; i<=rlength; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i],rf_fib[i]);
  }
  fclose(broken);
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
		       double **y_cp,int npoints)
{
  int i;
  double dy;

  for (i = 1; i <=npoints; i++) quick_interp(r_cp,y_cp,r[i],y,i);
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

void resize_and_interp(double h,double ****c,double ***s,double ***y,double **r,
		       double **rf_fib,double **y_cp,double *r_cp,int npoints,
		       int last_npoints,struct arr_ns *ns)
{
  copy_2_arrays(*r,*y,r_cp,y_cp,last_npoints); // copy arrays r and y into r_cp and y_cp
  free_matrices(c,s,y,r,rf_fib, // resize all arrays to new npoints size
		last_npoints,ns);
  allocate_matrices(c,s,y,r,rf_fib,
		    npoints,ns);
  propagate_r(*r,h,npoints);
  interpolate_array(*r,*y,r_cp,y_cp,npoints); // interpolate old y values (now stored in y_cp)
  //                                             so that the new y array has an initial guess
  //                                             for solvde.
  

}

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_fib,int npoints,struct arr_ns *ns)
{
  *y = matrix(1,ns->nyj,1,npoints);
  *s = matrix(1,ns->nsi,1,ns->nsj);
  *c = f3tensor(1,ns->nci,1,ns->ncj,1,npoints+1);
  *r = vector(1,npoints);
  *rf_fib = vector(1,npoints);
  return;
}

void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_fib,int npoints,struct arr_ns *ns)
{

  free_f3tensor(*c,1,ns->nci,1,ns->ncj,1,npoints+1);
  free_matrix(*s,1,ns->nsi,1,ns->nsj);
  free_matrix(*y,1,ns->nyj,1,npoints);
  free_vector(*r,1,npoints);
  free_vector(*rf_fib,1,npoints);
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
		 int last_npoints)
{
  int i;

  for (i = 1; i <=last_npoints; i++) {
    y_cp[1][i] = y[1][i];
    y_cp[2][i] = y[2][i];
    r_cp[i] = r[i];
  }

  return;
}

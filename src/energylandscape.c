#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"



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

void write_failure(double *r, double **y,double *rf_,int rlength,char *f_err);
void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt);
void propagate_r(double *r, double h,int mpt);
void save_psi(FILE *psi,double *r, double **y,int mpt);
void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable);
void make_f_err(char *f_err,int f_err_size,struct params p);
void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_npoints);
void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int npoints);
void quick_interp(double *xcp,double **ycp,double x,double **y,int i);
void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],double *R, 
			double *dEdR,double *dEdRlast,double *eta,
			double *dEdeta, double *dEdetalast,double *delta,
			double *dEddelta,double *dEddeltalast);

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns);
void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns);


void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt)
{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess!
    r[k] = (k-1)*h;
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

void save_psi(FILE *psi,double *r, double **y,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(psi,"%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i]);
  }
  printf("psi(R) = %1.2e\n",y[1][mpt]);
  return;
}

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable)
{
  fprintf(energy,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,E,
	  derivative,observable);
  return;
}

void make_f_err(char *f_err,int f_err_size,struct params p)
{
  snprintf(f_err,f_err_size,"data/QROMB_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,
	   p.delta,p.gamma_s);
  return;
}


void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
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

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int npoints)
{
  int i;
  double dy;

  for (i = 1; i <=npoints; i++) quick_interp(r_cp,y_cp,r[i],y,i);
  return;
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

void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],double *R, 
			double *dEdR,double *dEdRlast,double *eta,
			double *dEdeta, double *dEdetalast,double *delta,
			double *dEddelta,double *dEddeltalast)
// convenient way to make a variables (ending in var) that have the  //
// same addresses as whichever parameter that is specified by        //
// scan_what.
{

  if (strcmp(scan_what,"R")==0) {
    printf("R!\n");
    *var = R;
    *var0 = *R;
    *dEdvar = dEdR;
    *dEdvarlast = dEdRlast;
  }
  else if (strcmp(scan_what,"eta")==0) {
    printf("eta!\n");
    *var = eta;
    *var0 = *eta;
    *dEdvar = dEdeta;
    *dEdvarlast = dEdetalast;
  }
  else if (strcmp(scan_what,"delta")==0) {
    printf("delta!\n");
    *var = delta;
    *var0 = *delta;
    *dEdvar = dEddelta;
    *dEdvarlast = dEddeltalast;
  }
  else {
    printf("Need either R, eta, or delta as argv[n] input."
	   "Exiting to system.\n");
    exit(1);
  }
  return;
}


void write_failure(double *r, double **y,double *rf_,int rlength,char *f_err)
{
  int i;
  FILE *broken;
  printf("failed to integrate at npoints = %d, with R = %e.\n",rlength,r[rlength]);
  printf("saving psi(r) shape, and exiting to system.\n");
  broken = fopen(f_err,"w");
  for (i = 1; i<=rlength; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i],rf_[i]);
  }
  fclose(broken);
  exit(1);
  return;
}

void resize_all_arrays(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns)
{
  int last_npoints;

  last_npoints = (npoints-1)/2+1;
  free_matrices(c,s,y,r,rf_,integrand1,integrand2,last_npoints,ns);
  allocate_matrices(c,s,y,r,rf_,integrand1,integrand2,npoints,ns);


  return;
}

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns)
{
  *y = matrix(1,ns->nyj,1,npoints);
  *s = matrix(1,ns->nsi,1,ns->nsj);
  *c = f3tensor(1,ns->nci,1,ns->ncj,1,npoints+1);
  *r = vector(1,npoints);
  *rf_ = vector(1,npoints);
  *integrand1 = vector(1,npoints);
  *integrand2 = vector(1,npoints);
  return;
}

void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns)
{

  free_f3tensor(*c,1,ns->nci,1,ns->ncj,1,npoints+1);
  free_matrix(*s,1,ns->nsi,1,ns->nsj);
  free_matrix(*y,1,ns->nyj,1,npoints);
  free_vector(*r,1,npoints);
  free_vector(*rf_,1,npoints);
  free_vector(*integrand1,1,npoints);
  free_vector(*integrand2,1,npoints);
  return;
}


void scanE(struct params p,FILE *energy,FILE *psi,double conv,
	   int itmax,int mpt,int num_x,char scan_what[])
// The energy of the system E(R,L,eta). This function  //
// generates data of E vs scan_what[] (either "delta", //
// "R", or "eta"), while holding the other two values  //
// constant. The E vs scan_what[] data is saved in     //
// energy file. If a minimum (or multiple minima) are  //
// found in the energy landscape, psi(r) is saved to   //
// the psi file. //
{
  int npoints = mpt;
  int last_npoints = mpt;
  double *var,var0;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double dEdR, dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar, *dEdvarlast;
  double Emin = 1e100;
  int f_err_size = 200;
  char f_err[f_err_size];
  bool successful_calc = false;
  struct arr_ns ns;
  assign_ns(&ns);
  int max_size = (mpt-1)*4+1;
  double dx;
  int count_x;


  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);


  setup_var_pointers(&var,&var0,&dEdvar,&dEdvarlast,scan_what,&p.R, 
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,
		     &p.delta,&dEddelta,&dEddeltalast);
  dx = (p.upperbound_x-var0)/num_x;
  printf("dx = %lf\n",dx);
  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values
  
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  initialSlope = M_PI/(4.0*p.R);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  count_x = 0;

  while (count_x <= num_x) {

    // for each value of x (i.e *var) in E vs x

    do {


      h = p.R/(npoints-1);

      // only executes if the first calculation for the current var value is unsuccessful
      if (npoints != last_npoints) {
	printf("interpolating at %s = %e...\n",scan_what,*var);
	copy_arrays(r,y,r_cp,y_cp,last_npoints);
	free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		      last_npoints,&ns);
	allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
			  npoints,&ns);
	propagate_r(r,h,npoints);
	interpolate_array(r,y,r_cp,y_cp,npoints);
	printf("done interpolating.\n");

      }

      // if last loop was successful with last_xpoint sized arrays, then just
      // use last y values as initial guess
      else propagate_r(r,h,npoints);
      
      solvde(itmax,conv,slowc,scalv,&ns,npoints,y,r,c,s,&p,h); // relax to compute psi, psi' curves,

      // calculate energy, derivatives (see energy.c for code)
      successful_calc = energy_stuff(&E,&dEdR,&dEdeta,&dEddelta,&p,r,y,
				     rf_,integrand1,integrand2,npoints);

      npoints = (npoints-1)*2+1;

    }
    while (!successful_calc && npoints <= max_size);

    // since multiply npoints at the end of every loop, need to undo
    // one multiplication after the loop is finished.
    npoints = (npoints-1)/2+1;
    last_npoints = npoints;

    if (!successful_calc) { // writes the psi(r), r*f, etc, and exits
      make_f_err(f_err,f_err_size,p);
      write_failure(r,y,rf_,npoints,f_err);
    }

    // save var,E, and surface twist
    saveEnergy(energy,*var,E,*dEdvar,y[1][npoints]);
      


    // if a minimum has been found, save the psi(r) curve
    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0
	&& *dEdvarlast < 0 && E <= Emin) {
      save_psi(psi,r,y,npoints);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E+0.5);
      Emin = E;
    }

    count_x += 1;
    *var = var0+dx*count_x;
    dEdRlast = dEdR;
    dEdetalast = dEdeta;
    dEddeltalast = dEddelta;
    
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);

  return;
}



void scan2dE(struct params p,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,int num_x, int num_y,
	     char scan_what_x[],char scan_what_y[])
// The energy of the system E(R,L,eta). This function    //
// generates data of E vs scan_what_x[] vs scan_what_y[] //
// (either "delta","R", or "eta"), while holding the     //
// other value constant. The E vs x vs y data is saved   //
// in energy file. If a minimum (or multiple minima) are //
// found in the energy landscape, psi(r) are saved to    //
// the psi file. //
{
  int npoints = mpt;
  int last_npoints = mpt;
  double *var_x,var_x0;
  double *var_y,var_y0;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double dEdR, dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar_x, *dEdvar_xlast;
  double *dEdvar_y, *dEdvar_ylast;
  double Emin = 1e100;
  int f_err_size = 200;
  char f_err[f_err_size];
  bool successful_calc = false;
  struct arr_ns ns;
  assign_ns(&ns);
  int max_size = (mpt-1)*10+1;
  double dx,dy;
  int count_y,count_x;

  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  

  setup_var_pointers(&var_x,&var_x0,&dEdvar_x,&dEdvar_xlast,scan_what_x,&p.R, 
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,&p.delta,
		     &dEddelta,&dEddeltalast);
  dx = (p.upperbound_x-var_x0)/(num_x-1);

  setup_var_pointers(&var_y,&var_y0,&dEdvar_y,&dEdvar_ylast,scan_what_y,&p.R,
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,&p.delta,
		     &dEddelta,&dEddeltalast);
  dy = (p.upperbound_y-var_y0)/(num_y-1);

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  printf("dx = %lf, dy = %lf\n",dx,dy);
  
  initialSlope = M_PI/(4.0*p.R);

  count_y = 0;

  while (count_y < num_y) {
    npoints = mpt;
    last_npoints = mpt;
    count_x = 0;
    *var_x = var_x0;
    h = p.R/(npoints-1);
    // initial guess for the functional form of psi(r) and psi'(r)
    printf("initial slope for guess = %lf\n", initialSlope);
    linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

    
    while (count_x < num_x) {

      do {
	h = p.R/(npoints-1);
	
	// only executes if the first calculation for the current var value is unsuccessful
	if (npoints != last_npoints) {
	  printf("interpolating at %s = %e, %s = %e...\n",scan_what_x,*var_x,
		 scan_what_y,*var_y);
	  copy_arrays(r,y,r_cp,y_cp,last_npoints);
	  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
			last_npoints,&ns);
	  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
			    npoints,&ns);
	  propagate_r(r,h,npoints);
	  interpolate_array(r,y,r_cp,y_cp,npoints);
	  printf("done interpolating.\n");
	  
	}

	// if last loop was successful with last_xpoint sized arrays, then just
	// use last y values as initial guess
	else propagate_r(r,h,npoints);
	
	solvde(itmax,conv,slowc,scalv,&ns,npoints,y,r,c,s,&p,h); // relax to compute psi,
	//                                  psi' curves, 

	if (count_x == 0) initialSlope = 0.1*y[2][npoints];
	
	// calculate energy, derivatives (see energy.c for code)
	successful_calc = energy_stuff(&E,&dEdR,&dEdeta,&dEddelta,&p,r,y,
				       rf_,integrand1,integrand2,npoints);
      
      
	npoints = (npoints-1)*2+1;
	
      }
      while (!successful_calc && npoints <= max_size);
      // since multiply npoints at the end of every loop, need to undo
      // one multiplication after the loop is finished.
      npoints = (npoints-1)/2+1;
      last_npoints = npoints;

      if (!successful_calc) { // writes the psi(r), r*f, etc, and exits
	make_f_err(f_err,f_err_size,p);
	write_failure(r,y,rf_,npoints,f_err);
      }

      
      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",*dEdvar_x);
      fprintf(deriv_energy_y,"%.8e\t",*dEdvar_y);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);

      if (*var_x != var_x0 && *dEdvar_x*(*dEdvar_xlast) <= 0
	  && *dEdvar_xlast < 0 && *var_y != var_y0
	  && *dEdvar_y*(*dEdvar_ylast) <= 0 && *dEdvar_ylast <0) {
	save_psi(psi,r,y,npoints);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E+0.5);
	Emin = E;
      }

      dEdRlast = dEdR;
      dEdetalast = dEdeta;
      dEddeltalast = dEddelta;
      count_x += 1;
      *var_x = var_x0+count_x*dx;
    }
    fprintf(energy,"\n");
    fprintf(deriv_energy_x,"\n");
    fprintf(deriv_energy_y,"\n");
    fprintf(surfacetwist,"\n");
    count_y += 1;
    *var_y = var_y0+count_y*dy;
    printf("%s = %lf\n",scan_what_x,*var_x-dx);
    printf("%s = %lf\n",scan_what_y,*var_y-dy);
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  
  return;
}

void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *deriv_energy_x,FILE *deriv_energy_y,
	      FILE *surfacetwist,double conv,int itmax,
	      int mpt,double rateR, double rateeta,
	      double ratedelta)

{
  int npoints = mpt;
  int last_npoints = mpt;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double dEdR, dEdeta,dEddelta;
  int f_err_size = 200;
  char f_err[f_err_size];
  bool successful_calc = false;
  struct arr_ns ns;
  assign_ns(&ns);
  int max_size = (mpt-1)*4+1;
  double Rdot,etadot,deltadot;

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  
  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values
  
  initialSlope = M_PI/(4.0*p.R);
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  printf("initial slope for guess = %lf\n", initialSlope);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  dEdR = dEdeta = dEddelta = 1e100;

  // using classical gradient descent, try to find minimum.

  while (fabs(dEdR) > 1e-8 && fabs(dEdeta) > 1e-8
	 && fabs(dEddelta) > 1e-8) {

    do {
      h = p.R/(npoints-1);
      
      // only executes if the first calculation for the current variable values is unsuccessful
      if (npoints != last_npoints) {
	printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	       p.R,p.eta,p.delta);
	copy_arrays(r,y,r_cp,y_cp,last_npoints);
	free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		      last_npoints,&ns);
	allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
			  npoints,&ns);
	propagate_r(r,h,npoints);
	interpolate_array(r,y,r_cp,y_cp,npoints);
	printf("done interpolating.\n");
	
      }
      
      // if last loop was successful with last_xpoint sized arrays, then just
      // use last y values as initial guess
      else propagate_r(r,h,npoints);
	
      solvde(itmax,conv,slowc,scalv,&ns,npoints,y,r,c,s,&p,h); // relax to compute psi,
      //                                  psi' curves, 
      
      // calculate energy, derivatives (see energy.c for code)
      successful_calc = energy_stuff(&E,&dEdR,&dEdeta,&dEddelta,&p,r,y,
				     rf_,integrand1,integrand2,npoints);
      
      
      npoints = (npoints-1)*2+1;
      
    }
    while (!successful_calc && npoints <= max_size);
    // since multiply npoints at the end of every loop, need to undo
    // one multiplication after the loop is finished.
    npoints = (npoints-1)/2+1;
    last_npoints = npoints;
    
    if (!successful_calc) { // writes the psi(r), r*f, etc, and exits
      make_f_err(f_err,f_err_size,p);
      write_failure(r,y,rf_,npoints,f_err);
    }

    fprintf(energy,"%.8e\t",E);
    //    fprintf(deriv_energy_x,"%.8e\t",*dEdvar_x);
    //fprintf(deriv_energy_y,"%.8e\t",*dEdvar_y);
    fprintf(surfacetwist,"%.8e\t",y[1][mpt]);


    Rdot = -rateR*dEdR;
    etadot = -rateeta*dEdeta;
    deltadot = -ratedelta*dEddelta;

    p.R = p.R+Rdot;
    p.eta = p.eta+etadot;
    p.delta = p.delta+deltadot;

  }

  save_psi(psi,r,y,npoints);
  printf("SAVED!\n");
  printf("E_min-E_chol = %1.2e\n",E+0.5);

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  
  return;
}


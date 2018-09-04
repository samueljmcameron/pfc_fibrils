#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"


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

void make_f_err(char *f_err,int f_err_size,double K33,double k24,
		double Lambda,double d0,double omega,double R,
		double eta,double delta,double gamma_s)
{
  snprintf(f_err,f_err_size,"data/QROMB_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   K33,k24,Lambda,d0,omega,R,eta,delta,gamma_s);
  return;
}

void scanE(double K33,double k24,double Lambda,double d0,
	   double omega,double R,double eta,double delta,
	   double gamma_s,FILE *energy,FILE *psi,
	   double conv,int itmax,int mpt,
	   double upperbound,char scan_what[])
// The energy of the system E(R,L,eta). This function  //
// generates data of E vs scan_what[] (either "delta", //
// "R", or "eta"), while holding the other two values  //
// constant. The E vs scan_what[] data is saved in     //
// energy file. If a minimum (or multiple minima) are  //
// found in the energy landscape, psi(r) is saved to   //
// the psi file. //
{
  int xpoints = mpt;
  int last_xpoints = mpt;
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
  int ne = 2; // number of equations in ODE (2)
  int nb = 1; // number of BCs at first boundary
  int nsi = ne;
  int nsj = 2*ne+1;
  int nyj = ne;
  int nci = ne;
  int ncj = ne-nb+1;
  int max_size = (mpt-1)*4+1;

  y = matrix(1,nyj,1,mpt);
  y_cp = matrix(1,nyj,1,max_size);
  s = matrix(1,nsi,1,nsj);
  c = f3tensor(1,nci,1,ncj,1,mpt+1);
  r = vector(1,mpt);
  r_cp = vector(1,max_size);
  rf_ = vector(1,mpt);
  integrand1 = vector(1,mpt);
  integrand2 = vector(1,mpt);


  setup_var_pointers(&var,&var0,&dEdvar,&dEdvarlast,scan_what,&R, 
		     &dEdR,&dEdRlast,&eta,&dEdeta,&dEdetalast,
		     &delta,&dEddelta,&dEddeltalast);

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values
  

  while (*var <= upperbound) {

    // for each value of x (i.e *var) in E vs x

    do {

      h = R/(xpoints-1);

      // only executes if the first calculation for the current var value is unsuccessful
      if (xpoints != last_xpoints) {
	printf("interpolating...\n");
	copy_arrays(r,y,r_cp,y_cp,last_xpoints);
	resize_all_arrays(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
      			  xpoints,nci,ncj,nsi,nsj,nyj);
	propagate_r(r,h,xpoints);
	interpolate_array(r,y,r_cp,y_cp,xpoints);
	printf("done interpolating.\n");

      }
      // only executes the first loop, as after that any unsuccessful calculations exit
      // to the system.
      if (!successful_calc) { 
	initialSlope = M_PI/(2.01*R);
	linearGuess(r,y,initialSlope,h,xpoints); //linear initial guess 
      }
      // if last loop was successful with last_xpoint sized arrays, then just
      // use last y values as initial guess
      else propagate_r(r,h,xpoints);
      
      solvde(itmax,conv,slowc,scalv,ne,nb,xpoints,y,r,c,s,K33,k24,
	     Lambda,d0,eta,delta,h); // relax to compute psi, psi' curves,

      // calculate energy, derivatives (see energy.c for code)
      successful_calc = energy_stuff(&E,&dEdR,&dEdeta,&dEddelta,K33,k24,
				     Lambda,d0,omega,R,eta,delta,gamma_s,
				     r,y,rf_,integrand1,integrand2,xpoints);

      xpoints = (xpoints-1)*2+1;

    }
    while (!successful_calc && xpoints <= max_size);

    // since multiply xpoints at the end of every loop, need to undo
    // one multiplication after the loop is finished.
    xpoints = (xpoints-1)/2+1;
    last_xpoints = xpoints;

    if (!successful_calc) { // writes the psi(r), r*f, etc, and exits
      make_f_err(f_err,f_err_size,K33,k24,Lambda,d0,omega,R,
		 eta,delta,gamma_s);
      write_failure(r,y,rf_,xpoints,f_err);
    }

    // save var,E, and surface twist
    saveEnergy(energy,*var,E,*dEdvar,y[1][xpoints]);
      


    // if a minimum has been found, save the psi(r) curve
    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0
	&& *dEdvarlast < 0 && E <= Emin) {
      save_psi(psi,r,y,xpoints);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E+0.5);
      Emin = E;
    }

    *var += 0.001;
    dEdRlast = dEdR;
    dEdetalast = dEdeta;
    dEddeltalast = dEddelta;
    
  }

  free_f3tensor(c,1,nci,1,ncj,1,xpoints+1);
  free_matrix(s,1,nsi,1,nsj);
  free_matrix(y,1,nyj,1,xpoints);
  free_matrix(y_cp,1,nyj,1,max_size);
  free_vector(r,1,xpoints);
  free_vector(r_cp,1,max_size);
  free_vector(rf_,1,xpoints);
  free_vector(integrand1,1,xpoints);
  free_vector(integrand2,1,xpoints);

  return;
}

void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_xpoints)
{
  int i;

  for (i = 1; i <=last_xpoints; i++) {
    y_cp[1][i] = y[1][i];
    y_cp[2][i] = y[2][i];
    r_cp[i] = r[i];
  }

  return;
}

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int xpoints)
{
  int i;
  int last_xpoints = (xpoints-1)/2+1;
  double dy;

  for (i = 1; i <=xpoints; i++) {
    quick_interp(r_cp,y_cp[1],r[i],&y[1][i],last_xpoints);
    quick_interp(r_cp,y_cp[2],r[i],&y[2][i],last_xpoints);
  }
  return;
}

void quick_interp(double *xa, double *ya, double x, double *y,
		  int xpoints)
{
  int i=0;
  double diff=0.5*(xa[2]-xa[1]);
  int index;

  do i += 1;
  while (fabs(x-xa[i])>diff);
  if (i > xpoints) {printf("too many points!\n"); exit(1);}
  if (x==xa[i]) { *y = y[i]; }
  else if (x<xa[i]) i -= 1;
  double slope = (ya[i+1]-ya[i])/(2*diff);
  *y = slope*(x-xa[i+1])+ya[i+1];
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


void scan2dE(double *r,double **y,double ***c,double **s,
	     double K33,double k24,double Lambda,double d0,
	     double omega,double R,double eta,double delta,
	     double gamma_s,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,double upperbound_x,double upperbound_y,
	     char scan_what_x[],char scan_what_y[])
// The energy of the system E(R,L,eta). This function    //
// generates data of E vs scan_what_x[] vs scan_what_y[] //
// (either "delta","R", or "eta"), while holding the     //
// other value constant. The E vs x vs y data is saved   //
// in energy file. If a minimum (or multiple minima) are //
// found in the energy landscape, psi(r) are saved to    //
// the psi file. //
{
  int isitone;
  double *var_x,var_x0;
  double *var_y,var_y0;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double initialSlope;
  double rf_[mpt+1],integrand1[mpt+1],integrand2[mpt+1];
  double E;
  double dEdR,  dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar_x, *dEdvar_xlast;
  double *dEdvar_y, *dEdvar_ylast;
  double Emin = 1e100;
  int f_err_size = 200;
  char f_err[f_err_size];


  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.
  setup_var_pointers(&var_x,&var_x0,&dEdvar_x,&dEdvar_xlast,scan_what_x,&R, 
		     &dEdR,&dEdRlast,&eta,&dEdeta,&dEdetalast,&delta,
		     &dEddelta,&dEddeltalast);

  setup_var_pointers(&var_y,&var_y0,&dEdvar_y,&dEdvar_ylast,scan_what_y,&R, 
		     &dEdR,&dEdRlast,&eta,&dEdeta,&dEdetalast,&delta,
		     &dEddelta,&dEddeltalast);

  
  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values
  

  while (*var_x <= upperbound_x) {
    isitone = 1;
    *var_y = var_y0;
    while (*var_y <= upperbound_y) {

      h = R/(mpt-1);



      if (isitone == 1) {
	initialSlope = M_PI/(4.0*R);
	linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 
	isitone += 1;
      }
    
      else propagate_r(r,h,mpt); // if not first loop, previous 
      //                            result is initial guess

      make_f_err(f_err,f_err_size,K33,k24,Lambda,d0,omega,R,
		 eta,delta,gamma_s);

      solvde(itmax,conv,slowc,scalv,2,1,mpt,y,r,c,s,K33,k24,
	     Lambda,d0,eta,delta,h); // relax to compute psi,
      //                                  psi' curves, note the 2,1
      //                                  corresponds to two eqns,
      //                                   and 1 BC at the r = 0.

      //      if (*var_x > 0.002) {
      //	printf("%s = %lf\n",scan_what_x,*var_x);
      //	printf("%s = %lf\n",scan_what_y,*var_y);
      //	save_psi(psi,r,y,mpt);
      //  }
      // calculate energy, derivatives (see energy.c for code)
      energy_stuff(&E,&dEdR,&dEdeta,&dEddelta,K33,k24,Lambda,
		   d0,omega,R,eta,delta,gamma_s,r,y,rf_,
		   integrand1,integrand2,mpt);
      
      

      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",*dEdvar_x);
      fprintf(deriv_energy_y,"%.8e\t",*dEdvar_y);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);
      
      if (*var_x != var_x0 && *dEdvar_x*(*dEdvar_xlast) <= 0
	  && *dEdvar_xlast < 0 && *var_y != var_y0
	  && *dEdvar_y*(*dEdvar_ylast) <= 0 && *dEdvar_ylast <0) {
	save_psi(psi,r,y,mpt);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E+0.5);
	Emin = E;
      }

      *var_y += 0.001;
      dEdRlast = dEdR;
      dEdetalast = dEdeta;
      dEddeltalast = dEddelta;
    }
    fprintf(energy,"\n");
    fprintf(deriv_energy_x,"\n");
    fprintf(deriv_energy_y,"\n");
    fprintf(surfacetwist,"\n");
    printf("%s = %lf\n",scan_what_x,*var_x);
    printf("%s = %lf\n",scan_what_y,*var_y);
    *var_x += 0.001;
  }

  return;
}

void write_failure(double *r, double **y,double *rf_,int rlength,char *f_err)
{
  int i;
  FILE *broken;
  printf("failed to integrate at xpoints = %d, with R = %e.\n",rlength,r[rlength]);
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
		       double **integrand2,int xpoints,
		       int nci,int ncj,int nsi,int nsj,int nyj)
{
  int last_xpoints;

  last_xpoints = (xpoints-1)/2+1;
  free_f3tensor(*c,1,nci,1,ncj,1,last_xpoints+1);
  free_matrix(*s,1,nsi,1,nsj);
  free_matrix(*y,1,nyj,1,last_xpoints);
  free_vector(*r,1,last_xpoints);
  free_vector(*rf_,1,last_xpoints);
  free_vector(*integrand1,1,last_xpoints);
  free_vector(*integrand2,1,last_xpoints);
  
  *y = matrix(1,nyj,1,xpoints);
  *s = matrix(1,nsi,1,nsj);
  *c = f3tensor(1,nci,1,ncj,1,xpoints+1);
  *r = vector(1,xpoints);
  *rf_ = vector(1,xpoints);
  *integrand1 = vector(1,xpoints);
  *integrand2 = vector(1,xpoints);
  
  return;
}

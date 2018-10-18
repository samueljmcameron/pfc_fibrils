/* Functions for computing phase field crystal model of fibrils.


 */





#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"
#include <time.h>







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
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double *dEdx, *lastdEdx;
  double initialSlope;
  double E;
  double *dEdvar, *dEdvarlast;
  double Emin = 1e100;
  int max_size = (mpt-1)*4+1;
  double dx;
  int count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);

  setup_var_pointers(&var,&var0,&dEdvar,&dEdvarlast,scan_what,
		     &p,dEdx,lastdEdx);

  dx = (p.upperbound_x-var0)/num_x;
  printf("dx = %lf\n",dx);
  
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  initialSlope = M_PI/(4.0*p.R);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  count_x = 0;

  while (count_x <= num_x) {

    // for each value of x (i.e *var) in E vs x


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);

    last_npoints = npoints;

    // save var,E, and surface twist
    saveEnergy(energy,*var,E,*dEdvar,y[1][npoints]);
      


    // if a minimum has been found, save the psi(r) curve
    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0
	&& *dEdvarlast < 0 && E <= Emin) {
      save_psi(psi,r,y,npoints);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E);
      Emin = E;
    }

    count_x += 1;
    *var = var0+dx*count_x;
    arr_cp(lastdEdx,dEdx,3);
    
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);

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
  double *dEdx, *lastdEdx;
  double *dEdvar_x, *dEdvar_xlast;
  double *dEdvar_y, *dEdvar_ylast;
  double Emin = 1e100;
  int max_size = (mpt-1)*10+1;
  double dx,dy;
  int count_y,count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);

  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  

  setup_var_pointers(&var_x,&var_x0,&dEdvar_x,&dEdvar_xlast,scan_what_x,
		     &p,dEdx,lastdEdx);
  dx = (p.upperbound_x-var_x0)/(num_x-1);

  setup_var_pointers(&var_y,&var_y0,&dEdvar_y,&dEdvar_ylast,scan_what_y,
		     &p,dEdx,lastdEdx);
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

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		  y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);

      last_npoints = npoints;

      if (count_x == 0) initialSlope = 0.1*y[2][npoints];

      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",*dEdvar_x);
      fprintf(deriv_energy_y,"%.8e\t",*dEdvar_y);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);

      if (*var_x != var_x0 && *dEdvar_x*(*dEdvar_xlast) <= 0
	  && *dEdvar_xlast < 0 && *var_y != var_y0
	  && *dEdvar_y*(*dEdvar_ylast) <= 0 && *dEdvar_ylast <0) {
	save_psi(psi,r,y,npoints);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E);
	Emin = E;
      }

      arr_cp(lastdEdx,dEdx,3);
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
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);
  
  return;
}

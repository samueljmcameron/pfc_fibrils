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







void scanE(struct params p,double *x,FILE *energy,FILE *psi,
	   double conv,int itmax,int mpt,int num_scan,int x_index)
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
  double x0=x[x_index];
  double h;
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_fib;
  double *dEdx, *lastdEdx;
  double initialSlope;
  double E;
  double Emin = 1e100;
  int max_size = (mpt-1)*4+1;
  double scan_dx;
  int count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);


  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_fib,npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);

  scan_dx = (p.upperbound_x-x0)/num_scan;
  printf("scan_dx = %lf\n",scan_dx);
  
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  initialSlope = M_PI/(4.0*p.R);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  count_x = 0;

  while (count_x <= num_scan) {

    // for each value of x (i.e *var) in E vs x


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_fib,y_cp,r_cp,hessian,conv,itmax,
		&npoints,last_npoints,&ns,max_size);

    last_npoints = npoints;

    // save var,E, and surface twist
    saveEnergy(energy,x[x_index],E,dEdx[x_index],y[1][npoints]);
      


    // if a minimum has been found, save the psi(r) curve
    if (x[x_index] != x0 && dEdx[x_index]*lastdEdx[x_index] <= 0
	&& lastdEdx[x_index] < 0 && E <= Emin) {
      save_psi(psi,r,y,npoints);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E);
      Emin = E;
    }

    count_x += 1;
    x[x_index] = x0+scan_dx*count_x;
    arr_cp(lastdEdx,dEdx,3);
    
  }

  free_matrices(&c,&s,&y,&r,&rf_fib,npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);

  return;
}



void scan2dE(struct params p,double *x,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,int num_scanx, int num_scany,
	     int x_index, int y_index)
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
  double x0 = x[x_index];
  double y0 = x[y_index];
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_fib;
  double initialSlope;
  double E;
  double *dEdx, *lastdEdx;
  double Emin = 1e100;
  int max_size = (mpt-1)*10+1;
  double scan_dx,scan_dy;
  int count_y,count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);

  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  


  scan_dx = (p.upperbound_x-x0)/(num_scanx-1);

  scan_dy = (p.upperbound_y-y0)/(num_scany-1);

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  printf("scan_dx = %lf, scan_dy = %lf\n",scan_dx,scan_dy);
  
  initialSlope = M_PI/(4.0*x[1]);

  count_y = 0;

  while (count_y < num_scany) {
    npoints = mpt;
    last_npoints = mpt;
    count_x = 0;
    x[x_index] = x0;
    h = x[1]/(npoints-1);
    // initial guess for the functional form of psi(r) and psi'(r)
    printf("initial slope for guess = %lf\n", initialSlope);
    linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

    
    while (count_x < num_scanx) {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_fib,
		  y_cp,r_cp,hessian,conv,itmax,&npoints,
		  last_npoints,&ns,max_size,false);

      last_npoints = npoints;

      if (count_x == 0) initialSlope = 0.1*y[2][npoints];

      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",dEdx[x_index]);
      fprintf(deriv_energy_y,"%.8e\t",dEdx[y_index]);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);

      if (x[x_index] != x0 && dEdx[x_index]*(lastdEdx[x_index]) <= 0
	  && lastdEdx[x_index] < 0 && x[y_index] != y0
	  && dEdx[y_index]*(lastdEdx[y_index]) <= 0 && lastdEdx[y_index] <0) {
	save_psi(psi,r,y,npoints);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E);
	Emin = E;
      }

      arr_cp(lastdEdx,dEdx,3);
      count_x += 1;
      x[x_index] = x0+count_x*scan_dx;
    }
    fprintf(energy,"\n");
    fprintf(deriv_energy_x,"\n");
    fprintf(deriv_energy_y,"\n");
    fprintf(surfacetwist,"\n");
    count_y += 1;
    x[y_index] = y0+count_y*scan_dy;
    printf("x[%d] = %lf\n",x_index,x[x_index]-scan_dx);
    printf("x[%d] = %lf\n",y_index,x[y_index]-scan_dy);
  }

  free_matrices(&c,&s,&y,&r,&rf_,npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);
  
  return;
}

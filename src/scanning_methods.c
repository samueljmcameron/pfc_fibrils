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







void scanE(struct params p,double *x,FILE *energy,FILE *psi,const double conv,
	   const int itmax,int mpt,const int max_mpt,const int num_scan,
	   const int x_index,const int x_size)
/*==============================================================================

  Purpose: Scan across one out of three parameters in x = (R,eta,delta)' to map
  out E(x). The result will be saved data E vs x[x_index]

  ------------------------------------------------------------------------------

  Parameters:
n
  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector x = (R,eta,delta)' has the parameters which can be varied to
  get to equilibrium E(x).

  energy -- file where the E(x) vs x[x_index] and dEdx[x_index] vs x[x_index]
  will be saved.

  psi -- file where psi(r), and psi'(r) at the minimum E(x) in the 1d scan (so
  NOT the equilibrium psi(r)) are saved


  surfacetwist -- file where the surface twist vs x[x_index] is saved.

  conv -- convergence criterion for both solving ODE for psi(r) (once psi
  changes less than conv at each r) and when dEdx is small enough to say E(x)
  is at a minimum.

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  mpt -- starting grid size for r, psi(r) (interpolation is allowed, adding
  grid points. Note that mpt is of the form 2^n+1, where n is an integer.

  max_mpt -- maximum allowed grid size for r, psi(r) (still of the form 
  2^n+1)
  
  num_scan -- number of points in the x grid of the E(x) map.

  x_index -- value x[x_index] which will be scanned across.

  x_size -- size of the vector x. Since x = (R,eta,delta), this is 3.

  ------------------------------------------------------------------------------

  Returns: Does not return any values, but does save results into the files as
  described in the parameters section.

  ============================================================================*/
{

  void saveEnergy(FILE *energy, double R, double E, double derivative,
		  double observable);

  double **y;         // y[1][1..mpt] is psi(r), y[2][1..mpt] is psi'(r)
  double **y_cp;      // y_cp is copy of y
  double *r;          // r is radial distance from fibril centre
  double *r_cp;       // r_cp is a copy of r
  double *rf_fib;     // rf_fib is r*free energy density (integrand in model)
  double h;           // step spacing between points in r
  
  double **s,***c;     // dummy arrays for solvde func (ODE for psi(r))
  double initialSlope; // slope guess for very first psi(r) form
  
  double E;                       // E(x) energy of fibrils (cost function)
  double Emin = 1e100;            // E(x) at previous value of x
  double *dEdx, *lastdEdx;        // gradient vectors for E
  double *direction;              // direction of descent
  double *lastx;                  // will store a copy of x at its previous value
  double *x_scale;                // (R/p.Rscale,eta/p.etascale,delta/p.deltascale)'
  double *hessian;                // flattened hessian matrix of E
  double dx;                      // spacing used to calculate derivatives
  double cushion = 10.0;
  double *E_p,*E_m,*E_pij,*E_mij; // dummy matrices passed to derivative calcs

  struct arr_ns ns;         // used to set array sizes (for e.g. y, c, s)

  double x0 = x[x_index];

  double scan_dx;
  int count_x;


  assign_ns(&ns);

  scan_dx = (p.upperbound_x-x0)/num_scan;
  printf("scan_dx = %lf\n",scan_dx);
  

  // malloc the relevant arrays 
  allocate_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);

  // malloc the vectors of size x_size or x_size*x_size
  allocate_vectors(x_size,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
		   &E_p,&E_m,&E_pij,&E_mij);

  x_scale = vector(1,3);
  
  initialSlope = M_PI/(4.0*x[1]);
  h = x[1]/(mpt-1);

  linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 

  copy_2_arrays(r_cp,y_cp,r,y,mpt);

  array_constant(1e100,dEdx,x_size);
  arr_cp(lastdEdx,dEdx,x_size);


  count_x = 0;

  while (count_x <= num_scan) {

    x_scale[1] = x[1]/p.Rscale;
    x_scale[2] = x[2]/p.etascale;
    x_scale[3] = x[3]/p.deltascale;
    
    // for each value of x (i.e *var) in E vs x

    E = F_calc(&p,x_scale,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_mpt);

    dx = compute_dx(E,1e-8,cushion);

    derivatives_fd(dEdx,E,&p,x_scale,c,s,r,y,rf_fib,r_cp,y_cp,conv,dx,itmax,&mpt,&ns,
		   max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);

    // save var,E, and surface twist
    saveEnergy(energy,x[x_index],E,dEdx[x_index],y[1][mpt]);
      


    // if a minimum has been found, save the psi(r) curve
    if (positive_definite(hessian,x_size) && !non_zero_array(dEdx,conv,x_size)) {
      save_psi(psi,r,y,rf_fib,mpt);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E);
      Emin = E;
    }

    count_x += 1;
    x[x_index] = x0+scan_dx*count_x;
    arr_cp(lastdEdx,dEdx,3);
    
  }


  free_vectors(x_size,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
	       &E_p,&E_m,&E_pij,&E_mij);

  free_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);

  free_vector(x_scale,1,3);

  return;
}



void scan2dE(struct params p,double *x,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,FILE *surfacetwist,
	     const double conv,const int itmax,int mpt,const int max_mpt,
	     const int num_scanx,const int num_scany,const int x_index,
	     const int y_index, const int x_size)
/*==============================================================================

  Purpose: Scan across two out of three parameters in x = (R,eta,delta)' to map
  out E(x). The result will be saved data of the map of E vs x[x_index] vs
  x[y_index].

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector x = (R,eta,delta)' has the parameters which can be varied to
  get to equilibrium E(x).

  energy -- file where the map of E(x) vs x vs y is saved.

  psi -- file where psi(r), and psi'(r) at the minimum E(x) in the 2d scan (so
  NOT the equilibrium psi(r)) are saved

  deriv_energy_x, deriv_energy_y -- file where the maps of dEdx vs x vs y and
  dEdy vs x vs y are saved.

  surfacetwist -- file where the map of surface twist vs x vs y is saved.

  conv -- convergence criterion for both solving ODE for psi(r) (once psi
  changes less than conv at each r) and when dEdx is small enough to say E(x)
  is at a minimum.

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  mpt -- starting grid size for r, psi(r) (interpolation is allowed, adding
  grid points. Note that mpt is of the form 2^n+1, where n is an integer.

  max_mpt -- maximum allowed grid size for r, psi(r) (still of the form 
  2^n+1)
  
  num_scanx (num_scany) -- number of points in the x (y) grid of the E(x) map.

  x_index (y_index) -- value x[x_index] (x[y_index]) which will be scanned
  across.

  x_size -- size of the vector x. Since x = (R,eta,delta), this is 3.

  ------------------------------------------------------------------------------

  Returns: Does not return any values, but does save results into the files as
  described in the parameters section.

  ============================================================================*/
{



  double **y;         // y[1][1..mpt] is psi(r), y[2][1..mpt] is psi'(r)
  double **y_cp;      // y_cp is copy of y
  double *r;          // r is radial distance from fibril centre
  double *r_cp;       // r_cp is a copy of r
  double *rf_fib;     // rf_fib is r*free energy density (integrand in model)
  double h;           // step spacing between points in r
  
  double **s,***c;     // dummy arrays for solvde func (ODE for psi(r))
  double initialSlope; // slope guess for first psi(r) form when count_x == 0.
  
  double E;                       // E(x) energy of fibrils (cost function)
  double Emin = 1e100;            // E(x) at previous value of x
  double *dEdx, *lastdEdx;        // gradient vectors for E
  double *direction;              // direction of descent
  double *lastx;                  // will store a copy of x at its previous value
  double *x_scale;                // (R/p.Rscale,eta/p.etascale,delta/p.deltascale)'
  double *hessian;                // flattened hessian matrix of E
  double dx;                      // spacing used to calculate derivatives
  double cushion = 10;            //
  double *E_p,*E_m,*E_pij,*E_mij; // dummy matrices passed to derivative calcs

  struct arr_ns ns;         // used to set array sizes (for e.g. y, c, s)

  double x0 = x[x_index];
  double y0 = x[y_index];

  double scan_dx,scan_dy;
  int count_y,count_x;


  assign_ns(&ns);


  scan_dx = (p.upperbound_x-x0)/(num_scanx-1);
  scan_dy = (p.upperbound_y-y0)/(num_scany-1);
  printf("scan_dx = %lf, scan_dy = %lf\n",scan_dx,scan_dy);


  // malloc the relevant arrays 
  allocate_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);

  // malloc the vectors of size x_size or x_size*x_size
  allocate_vectors(x_size,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
		   &E_p,&E_m,&E_pij,&E_mij);

  x_scale = vector(1,3);
  
  array_constant(1e100,dEdx,x_size);
  arr_cp(lastdEdx,dEdx,x_size);

  initialSlope = M_PI/(4.0*x[1]);

  count_y = 0;

  while (count_y < num_scany) {

    count_x = 0;
    x[x_index] = x0;

    x_scale[1] = x[1]/p.Rscale;
    x_scale[2] = x[2]/p.etascale;
    x_scale[3] = x[3]/p.deltascale;

    printf("initial slope for guess = %lf\n", initialSlope);

    h = x[1]/(mpt-1);
    
    linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 
    
    copy_2_arrays(r_cp,y_cp,r,y,mpt);

    
    while (count_x < num_scanx) {

      E = F_calc(&p,x_scale,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_mpt);

      dx = compute_dx(E,1e-8,cushion);

      derivatives_fd(dEdx,E,&p,x_scale,c,s,r,y,rf_fib,r_cp,y_cp,conv,dx,itmax,&mpt,&ns,
		     max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);

      if (count_x == 0) initialSlope = 0.1*y[2][mpt];

      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",dEdx[x_index]);
      fprintf(deriv_energy_y,"%.8e\t",dEdx[y_index]);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);

      if (positive_definite(hessian,x_size) && !non_zero_array(dEdx,conv,x_size)) {
	save_psi(psi,r,y,rf_fib,mpt);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E);
	Emin = E;
      }

      arr_cp(lastdEdx,dEdx,x_size);
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


  free_vectors(x_size,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
	       &E_p,&E_m,&E_pij,&E_mij);

  free_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);

  free_vector(x_scale,1,3);

  
  return;
}

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable)
{
  fprintf(energy,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,E,
	  derivative,observable);
  return;
}

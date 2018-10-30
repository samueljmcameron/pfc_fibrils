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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


#define EFFECTIVE_ZERO 1e-14


void graddesc(struct params p,double *x,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,FILE *denergyddelta,
	      FILE *surfacetwist,FILE *energydensity,const double conv,
	      const int itmax,int mpt,const int max_mpt,double rate,
	      const int x_size0)

/*==============================================================================
  
  Purpose: This function minimizes the free energy functional E(x;psi(r)),
  where x = (R,eta,delta)', and psi(r) is the double-twist director field. It 
  does this by a combination of conjugate gradient descent and second order
  Newton-Raphson method. At each x, the form of psi(r) is determined by solving
  the ODE generated from setting the functional derivate dEdpsi = 0 with a 
  relaxation method (described in Chapter 17 of Numerical Recipes in C). The
  the gradient dEdx used in conjugate-gradient and Newton-Raphson is determined
  through finite-differences calculations, as dEdx cannot be found
  analytically.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the parameters to be optimized with descent, 
  where x = (R,eta,delta)'.

  energy -- file where the progression of E(x) is saved.

  psi -- file where the equilibrium (final) r, psi(r), and psi'(r) are saved

  denergydR, denergydeta, denergyddelta -- file where the progressions of x and
  dEdx are saved (each component in a different file)

  surfacetwist -- file where progression of surface twist psi(R) is saved.

  energydensity -- file where the equilibrium (final) r*f_fibril (energy
  integrand) is saved.

  conv -- convergence criterion for both solving ODE for psi(r) (once psi
  changes less than conv at each r) and when dEdx is small enough to say E(x)
  is at a minimum.

  itmax -- maximum number of iterations allowed when solving ODE for psi(r).

  mpt -- starting grid size for r, psi(r) (interpolation is allowed, adding
  grid points. Note that mpt is of the form 2^n+1, where n is an integer.

  max_mpt -- maximum allowed grid size for r, psi(r) (still of the form 
  2^n+1)
  
  rate -- maximum rate of descent considered.

  x_size0 -- size of the vector x. Since x = (R,eta,delta), this is 3.

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
  double initialSlope; // slope guess for very first psi(r) form
  
  int x_size = x_size0;           // x_size can change to 1 if delta goes to zero
  double E;                       // E(x) energy of fibrils (cost function)
  double lastE = 1e100;           // E(x) at previous value of x
  double *dEdx, *lastdEdx;        // gradient vectors for E
  double *direction;              // direction of descent
  double *lastx;                  // will store a copy of x at its previous value
  double *hessian;                // flattened hessian matrix of E
  double *E_p,*E_m,*E_pij,*E_mij; // dummy matrices passed to derivative calcs

  const double min_rate = 1e-10;  // minimum value of rate

  int count = 0;           // count how many times we descend before minimum reached

  clock_t begin = clock(); // timing function stuff
  clock_t end;

  int pos_def_in_a_row = 0; // count number of positive definite hessians in a row

  struct arr_ns ns;         // used to set array sizes (for e.g. y, c, s)

  assign_ns(&ns);           // set struct values


  // malloc the relevant arrays 
  allocate_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);

  // malloc the vectors of size x_size0 or x_size0*x_size0
  allocate_vectors(x_size0,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
		   &E_p,&E_m,&E_pij,&E_mij);
  
  initialSlope = M_PI/(4.0*x[1]);
  h = x[1]/(mpt-1);

  linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 

  copy_2_arrays(r_cp,y_cp,r,y,mpt);
  array_constant(1e100,dEdx,x_size0);
  arr_cp(lastdEdx,dEdx,x_size0);
  array_constant(0,direction,x_size0);



  // using classical gradient descent newton raphson hybrid, try to find minimum.
  // start with conjugate gradient descent until hessian seems to look positive definite

  while (pos_def_in_a_row < 20) {

    arr_cp(lastx,x,x_size0);

    //    if (x[2] <= EFFECTIVE_ZERO) x_size = 1;
    //    else x_size = x_size0;

    E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_mpt);

    if (count % 100 != 0 || count == 0) {

      derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,conv,itmax,&mpt,&ns,
		     max_mpt,x_size,hessian,false,E_p,E_m,E_pij,E_mij);
      

      set_direction(direction,dEdx,lastdEdx,x_size);

      armijo_backtracker(rate,E,dEdx,direction,&p,x,r,y,rf_fib,c,s,r_cp,y_cp,
			 conv,itmax,&mpt,&ns,max_mpt,min_rate,x_size);

    } else {

      derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,conv,itmax,&mpt,&ns,
		     max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);



      if (positive_definite(hessian,x_size)) {

	printf("Newton Raphson method worked! number of positive " 
	       "definite hessians in a row = %d!\n",
	       pos_def_in_a_row);

	pos_def_in_a_row += 1;

      } else {
	printf("Newton's Raphson method didn't work. Resetting "
	       "the number of positive definite hessians in "
	       "a row back to %d.\n",pos_def_in_a_row);

	pos_def_in_a_row = 0;
	
      }
      set_direction(direction,dEdx,lastdEdx,x_size);
      
      armijo_backtracker(rate,E,dEdx,direction,&p,x,r,y,rf_fib,c,s,r_cp,y_cp,
			 conv,itmax,&mpt,&ns,max_mpt,min_rate,x_size);



    }

    printf("dEdx[1] = %e\t",dEdx[1]);
    printf("dEdx[2] = %e\t",dEdx[2]);
    printf("dEdx[3] = %e\n",dEdx[3]);

    fprintf(energy,"%d\t%.8e\n",count,E);

    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastx[1],E,dEdx[1],
	    hessian[1]);

    fprintf(denergydeta,"%.8e\t%.8e\n",lastx[2],dEdx[2]);

    fprintf(denergyddelta,"%.8e\t%.8e\n",lastx[3],dEdx[3]);

    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastx[1],y[1][mpt],
	    y[2][mpt]);

    arr_cp(lastdEdx,dEdx,x_size0);
    count += 1;
    printf("count = %d\n",count);

    if (x[1] <= 0) {
      printf("R is being driven to negative values, R = %e!\n",x[1]);
      x[1] = 1e-6;
    }
      
  }

  printf("\n\n\n\n\n\n"
	 "moving to (strictly) Newton-Raphson method!"
	 "\n\n\n\n\n\n");


  while (non_zero_array(dEdx,conv,x_size)) {

    arr_cp(lastx,x,x_size0);

    //    if (x[3] <= EFFECTIVE_ZERO) x_size = 1;
    //    else x_size = x_size0;

    E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_mpt);

    derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,conv,itmax,&mpt,&ns,
		   max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);

    printf("hessian[1][1] = %e\n",hessian[1]);
    
    if (!positive_definite(hessian,x_size)) {      
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      
      set_direction(direction,dEdx,lastdEdx,x_size);
      
      armijo_backtracker(rate,E,dEdx,direction,&p,x,r,y,rf_fib,c,s,r_cp,y_cp,
			 conv,itmax,&mpt,&ns,max_mpt,min_rate,x_size);
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      
    } else {
      
      hessian_update_x(x,hessian,dEdx,x_size);

    }



    fprintf(energy,"%d\t%.8e\n",count,E);

    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastx[1],E,dEdx[1],
	    hessian[1]);

    fprintf(denergydeta,"%.8e\t%.8e\n",lastx[2],dEdx[2]);

    fprintf(denergyddelta,"%.8e\t%.8e\n",lastx[3],dEdx[3]);

    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastx[1],y[1][mpt],
	    y[2][mpt]);


    printf("dEdx[1] = %e\t",dEdx[1]);
    printf("dEdx[2] = %e\t",dEdx[2]);
    printf("dEdx[3] = %e\n",dEdx[3]);


    arr_cp(lastdEdx,dEdx,x_size0);
    lastE = E;

    count += 1;
    printf("count = %d\n",count);


      
  }

  E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_mpt);

  derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,conv,itmax,&mpt,&ns,
		 max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);

  positive_definite(hessian,x_size);

  end = clock();

  printf("calculation for %d iterations with Newton descent = %e\n",
	 count,(double) (end-begin)/CLOCKS_PER_SEC);


  printf("count = %d\n",count);
  save_psi(psi,r,y,mpt);
  save_energydensity(energydensity,r,rf_fib,mpt);
  printf("SAVED!\n");
  printf("E_min-E_chol = %1.2e\n",E);


  free_vectors(x_size0,&lastx,&dEdx,&lastdEdx,&direction,&hessian,
	       &E_p,&E_m,&E_pij,&E_mij);

  free_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,max_mpt);



  return;
}



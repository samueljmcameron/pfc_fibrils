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

#define HESS(i,j) hessian[(j)+(i-1)*(3)]


void graddesc(struct params p,double *x,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,FILE *denergyddelta,
	      FILE *surfacetwist,FILE *energydensity,const double convODE,
	      const double convMIN,const int itmax,int mpt,const int max_mpt,
	      double rate,const int x_size0)

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

  convODE -- convergence criterion for solving ODE for psi(r) (once psi
  changes less than convODE at each r).

  convMIN -- convergence criterion for when dEdx is small enough to say E(x)
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

  void print_x_dEdx(double *x, double *dEdx,double *hessian,double x_size);

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
  double dx;                      // spacing used when calculating finite differences
  double *E_p,*E_m,*E_pij,*E_mij; // dummy matrices passed to derivative calcs
  bool calc_hess = false;

  const double min_rate = 1e3*EFFECTIVE_ZERO;  // minimum value of rate

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

  calc_hess = true;
  while (non_zero_array(dEdx,convMIN,x_size)) {


    arr_cp(lastx,x,x_size0);

    x_size = 3;

    printf("start of loop,x[1] = %e\n",x[1]);
    E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,&mpt,&ns,max_mpt);

    printf("calculated E_calc, x = (%e,%e,%e)\n",x[1],x[2],x[3]);

    printf("E = %e\n",E);
    dx = compute_dx(E,convMIN);
    dx = (x[1]-dx>0) ? dx : x[1]-x[1]/2.0;

    printf("dx = %e\n",dx);


    derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,convODE,dx,itmax,
		   &mpt,&ns,max_mpt,x_size,hessian,calc_hess,E_p,E_m,E_pij,
		   E_mij);

    printf("passed derivative calc\n");
     
    if (fabs(x[3])<=convMIN) {
      x_size = 1;
    } else {
      x_size = 3;
    }
   
    if (calc_hess) {
      
      if (positive_definite(hessian,x_size)) {

	printf("Newton Raphson method worked! number of positive " 
	       "definite hessians in a row = %d!\n",
	       pos_def_in_a_row);

	//hessian_update_x(x,hessian,dEdx,x_size);

	pos_def_in_a_row += 1;
      
      } else {

	pos_def_in_a_row = 0;

	printf("Newton's Raphson method didn't work. Resetting "
	       "the number of positive definite hessians in "
	       "a row back to %d.\n",pos_def_in_a_row);


      }
    }

    set_direction(direction,dEdx,lastdEdx,x_size);    
	
    armijo_backtracker(rate,E,dEdx,direction,&p,x,r,y,rf_fib,c,s,r_cp,y_cp,
		       convODE,itmax,&mpt,&ns,max_mpt,min_rate,x_size);
	
    printf("direction[1] = %e\n",direction[1]);

    fprintf(energy,"%d\t%.8e\n",count,E);

    fprintf(denergydR,"%.8e\t%.8e\t%.8e\t%.8e\n",lastx[1],dEdx[1],
	    HESS(1,1),dx);

    fprintf(denergydeta,"%.8e\t%.8e\t%.8e\n",lastx[2],dEdx[2],HESS(2,2));

    fprintf(denergyddelta,"%.8e\t%.8e\t%.8e\n",lastx[3],dEdx[3],HESS(3,3));

    fprintf(surfacetwist,"%.8e\t%.8e\t%.8e\n",lastx[1],y[1][mpt],
	    y[2][mpt]);

    arr_cp(lastdEdx,dEdx,x_size0);
    count += 1;
    printf("count = %d, x_size = %d\n",count,x_size);
      
  }

  printf("\n\n\n\n\n\n"
	 "moving to (strictly) Newton-Raphson method!"
	 "\n\n\n\n\n\n");
  

  while (non_zero_array(dEdx,convMIN,x_size)) {



    x_size = 3;

    arr_cp(lastx,x,x_size0);

    E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,&mpt,&ns,max_mpt);

    dx = compute_dx(E,convMIN);

    derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,convODE,dx,itmax,
		   &mpt,&ns,max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);
    

    if (fabs(x[3])<=convMIN) {
      x_size = 1;
    } else {
      x_size = 3;
    }

    if (!positive_definite(hessian,x_size)) {      
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      
      set_direction(direction,dEdx,lastdEdx,x_size);
      
      armijo_backtracker(rate,E,dEdx,direction,&p,x,r,y,rf_fib,c,s,r_cp,y_cp,
      			 convODE,itmax,&mpt,&ns,max_mpt,min_rate,x_size);
      
      
    } else {
      
      hessian_update_x(x,hessian,dEdx,x_size);

    }



    fprintf(energy,"%d\t%.8e\n",count,E);

    fprintf(denergydR,"%.12e\t%.12e\t%.12e\t%.12e\n",lastx[1],dEdx[1],
	    HESS(1,1),dx);

    fprintf(denergydeta,"%.8e\t%.8e\t%.8e\n",lastx[2],dEdx[2],HESS(2,2));

    fprintf(denergyddelta,"%.8e\t%.8e\t%.8e\n",lastx[3],dEdx[3],HESS(3,3));

    fprintf(surfacetwist,"%.8e\t%.8e\t%.8e\n",lastx[1],y[1][mpt],
	    y[2][mpt]);


    arr_cp(lastdEdx,dEdx,x_size0);
    lastE = E;

    count += 1;
    printf("count = %d, x_size = %d\n",count,x_size);


      
  }

  arr_cp(lastx,x,x_size0);

  E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,convODE,itmax,&mpt,&ns,max_mpt);

  dx = compute_dx(E,convMIN);

  derivatives_fd(dEdx,E,&p,x,c,s,r,y,rf_fib,r_cp,y_cp,convODE,dx,itmax,
		 &mpt,&ns,max_mpt,x_size,hessian,true,E_p,E_m,E_pij,E_mij);


  fprintf(energy,"%d\t%.8e\n",count,E);

  fprintf(denergydR,"%.8e\t%.8e\t%.8e\t%.8e\n",lastx[1],dEdx[1],
	  HESS(1,1),dx);

  fprintf(denergydeta,"%.8e\t%.8e\t%.8e\n",lastx[2],dEdx[2],HESS(2,2));
  
  fprintf(denergyddelta,"%.8e\t%.8e\t%.8e\n",lastx[3],dEdx[3],HESS(3,3));
  
  fprintf(surfacetwist,"%.8e\t%.8e\t%.8e\n",lastx[1],y[1][mpt],
	  y[2][mpt]);

  print_x_dEdx(x,dEdx,hessian,x_size);

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


void print_x_dEdx(double *x, double *dEdx,double *hessian,double x_size)
{
  int i,j;

  //  for (i = 1; i <= x_size; i++) {
      //    printf("x[%d] = %e\tdEdx[%d] = %e\tHESS[%d][%d] = %e\n",
      //   i,x[i],i,dEdx[i],i,i,HESS(i,i));
  // }

  for (i = 1; i <= x_size; i++) {
    for (j = 1; j <= x_size; j++) {
      printf("%e\t",HESS(i,j));
    }
    printf("\n");
  }
  return;

}

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
	      const int itmax,int mpt,const int max_mpt,double rate0,
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
  
  rate0 -- maximum rate size in dR = -rate0*(descent direction).

  x_size0 -- size of the vector x. Since x = (R,eta,delta), this is 3.

  ------------------------------------------------------------------------------

  Returns: Does not return any values, but does save results into the files as
  described in the parameters section.

  ============================================================================*/
  
{


  double **y;         // y[1][1..mpt] is psi(r), y[2][1..mpt] is psi'(r)
  double **y_cp;      // y_cp is copy of y
  double *r,*r_cp;    // r is radial distance from fibril centre
  double *rf_fib;     // rf_fib is r*free energy density (integrand in model)
  double h;           // step spacing between points in r
  
  double **s,***c;     // dummy arrays for solvde func (ODE for psi(r))
  double initialSlope; // slope guess for very first psi(r) form
  int last_mpt = mpt;  // number of grid points in r and y

  int max_size = (mpt-1)*8+1; // maximum number of grid points for r and y
  
  int x_size = x_size0;           // x_size can change to 1 if delta goes to zero
  double *x_cp;                   // copy of x = (R,eta,delta)
  double E;                       // E(x) energy of fibrils (cost function)
  double *dEdx, *lastdEdx;        // gradient vectors for E
  double *hessian;                // flattened hessian matrix of E
  double *E_p,*E_m,*E_pij,*E_mij; // dummy matrices passed to derivative calcs

  double rate = rate0;            // current value of rate in dx = -rate*dEdx
  const double min_rate = 1e-10;  // minimum value of rate
  double betak;                   // parameter in conjugate gradient descent

  int count = 0;           // count how many times we descend before minimum reached

  clock_t begin = clock(); // timing function stuff
  clock_t end;

  int pos_def_in_a_row = 0; // count number of positive definite hessians in a row

  struct arr_ns ns;         // used to set array sizes (for e.g. y, c, s)

  assign_ns(&ns);           // set struct values


  // malloc the relevant arrays which may be resized
  allocate_matrices(&c,&s,&y,&r,&rf_fib,max_size,ns);

  // malloc the vectors of size x_size0 or x_size0*x_size0
  allocate_vectors(x_size0,&x,&dEdx,&lastdEdx,&direction,&hessian,
		   &x_cp,&E_p,&E_m,&E_pij,&E_mij);
  
  initialSlope = M_PI/(4.0*x[1]);
  h = x[1]/(mpt-1);
  // initial guess for the functional form of psi(r) and psi'(r)

  linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 
  linearGuess(r_cp,y_cp,initialSlope,h,mpt); //linear initial guess 

  // using classical gradient descent newton raphson hybrid, try to find minimum.


  array_constant(1e100,dEdx,x_size0);
  arr_cp(lastdEdx,dEdx,x_size0);
  array_constant(0,direction,x_size0);




  // start with conjugate gradient descent until hessian seems to look positive definite

  while (pos_def_in_a_row < 20 && fabs(dEdx[1])>conv) {

    if (x[2] <= EFFECTIVE_ZERO) x_size = 1;
    else x_size = x_size0;

    E = E_calc(&p,x,r,y,rf_fib,c,s,r_cp,y_cp,conv,itmax,&mpt,&ns,max_size);

    if (0==1) {//(count % 100 != 0 || count == 0)) {

      betak = polak_betak(dEdx,lastdEdx);

      set_direction(direction,dEdx,betak);

      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,mpt,&ns,last_mpt,max_size,
				min_rate,x_size);

      update_p(&p,rate,direction,dx,x_size);

    } else {


      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_fib,y_cp,r_cp,
		  hessian,conv,itmax,&mpt,last_mpt,&ns,max_size);

      if (positive_definite(hessian,x_size)) {

	hessian_update_p(&p,hessian,dEdx,dx,x_size);

	pos_def_in_a_row += 1;

	printf("Newton Raphson method worked! number of positive " 
	       "definite hessians in a row = %d!\n",
	       pos_def_in_a_row);

      } else {

	betak = polak_betak(dEdx,lastdEdx);
	
	set_direction(direction,dEdx,betak);

	rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				  conv,itmax,mpt,&ns,last_mpt,max_size,
				  min_rate,x_size);
	

	update_p(&p,rate,direction,dx,x_size);

	pos_def_in_a_row = 0;

	printf("Newton's Raphson method didn't work. Resetting "
	       "the number of positive definite hessians in "
	       "a row back to %d.\n",pos_def_in_a_row);

      }

    }
    
    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastR,E,dEdx[1],hessian[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastR,y[1][mpt],y[2][mpt]);

    last_mpt = mpt;
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    lastE = E;


    count += 1;
    printf("count = %d\n",count);

    if (p.R <= 0) {
      printf("R is being driven to negative values, R = %e!\n",p.R);
      p.R = 1e-6;
    }
      
  }

  printf("\n\n\n\n\n\n"
	 "moving to (strictly) Newton-Raphson method!"
	 "\n\n\n\n\n\n");


  while (non_zero_array(dEdx,conv)) {

    if (p.delta <= EFFECTIVE_ZERO) x_size = 1;
    else x_size = x_size0;

    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_fib,y_cp,r_cp,
		hessian,conv,itmax,&mpt,last_mpt,&ns,max_size);

    printf("hessian[1][1] = %e\n",hessian[1]);
    
    if (!positive_definite(hessian,x_size)) {      
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      betak = polak_betak(dEdx,lastdEdx);
      
      set_direction(direction,dEdx,betak);
      
      rate = armijo_backtracker(st,rate0,E,dEdx,direction,&p,r,y,r_cp,y_cp,
				conv,itmax,mpt,&ns,last_mpt,max_size,
				min_rate,x_size);
      
      update_p(&p,rate,direction,dx,x_size);
      
      printf("at count %d, energy got bigger by %e.\n",count,E-lastE);
      
    } else {
      
      hessian_update_p(&p,hessian,dEdx,dx,x_size);

    }


    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.14e\t%.14e\t%.14e\t%.14e\n",lastR,E,dEdx[1],hessian[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",lasteta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",lastdelta,dEdx[3]);
    fprintf(surfacetwist,"%.14e\t%.14e\t%.14e\n",lastR,y[1][mpt],y[2][mpt]);

    last_mpt = mpt;
    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    lastE = E;

    count += 1;
    printf("count = %d\n",count);


      
  }


  single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_fib,y_cp,r_cp,
	      hessian,conv,itmax,&mpt,last_mpt,
	      &ns,max_size,true);

  positive_definite(hessian,x_size);

  end = clock();

  printf("calculation for %d iterations with Newton descent = %e\n",
	 count,(double) (end-begin)/CLOCKS_PER_SEC);


  printf("count = %d\n",count);
  save_psi(psi,r,y,mpt);
  save_energydensity(energydensity,r,rf_fib,mpt);
  printf("SAVED!\n");
  printf("E_min-E_chol = %1.2e\n",E);

  free_matrices(&c,&s,&y,&r,&rf_fib,max_size,ns);



  return;
}



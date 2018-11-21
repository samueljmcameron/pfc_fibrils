#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"
#include <time.h>



void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x)
/*==============================================================================

  Purpose: This function writes a string to the char array f_err. f_err will be
  the file name of a file which stores e.g. psi(r) if the calculation of E(x)
  fails.

  ------------------------------------------------------------------------------

  Parameters:

  f_err -- The char array where the file name is stored.

  err_type -- A string that is inputted to the function, and written to the
  f_err file name. The string identifies the type of error that caused the
  calculation to fail. The two possible types of errors are that the ODE
  solver solvde_wrapper fails to find psi(r) - err_type = "SOLVDE", or that
  the calculation of E(x) fails as the integration cannot be done with a 
  sufficiently low error - err_type = "QROMB".

  f_err_size -- Approximate size of the f_err char array.

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  x -- This vector holds the variable parameters x = (R,eta,delta)'.

  ------------------------------------------------------------------------------

  Returns: Does not return a value, it just saves the file name to the f_err 
  char array.
  ============================================================================*/
  
{
  /*
  snprintf(f_err,f_err_size,"../../tmp_data/%s_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",err_type,
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,x[1],x[2],x[3],
	   p.gamma_s);
  */
  snprintf(f_err,f_err_size,"../../tmp_data/%s_psivsr.txt",err_type);
  return;
}


void save_psi(FILE *psi,double *r, double **y,double *rf_fib,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(psi,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i],rf_fib[i]);
  }
  printf("psi(R) = %1.2e\n",y[1][mpt]);
  return;
}

void save_energydensity(FILE *energydensity,double *r, double *rf_fib,
			int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(energydensity,"%10.8e\t%10.8e\n",r[i],rf_fib[i]);
  }
  return;
}

int x_index(char scan_what[])
{
  if (strcmp(scan_what,"R")==0) {
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

bool non_zero_array(double *dEdx,double conv,int x_size)
{
  int i;

  for (i = 1; i <= x_size; i++) {
    if (fabs(dEdx[i])>conv) return true;
  }
  return false;
}


void allocate_vectors(double x_size,double **lastx,double **dEdx,
		      double **lastdEdx,double **direction,double **hessian,
		      double **E_p,double **E_m,double **E_pij,double **E_mij)
/*==============================================================================
  ============================================================================*/

{

  // these vectors are meaningful
  *lastx = vector(1,x_size); // x = (R,eta,delta)
  *dEdx = vector(1,x_size); // dEdx = grad(E)
  *lastdEdx = vector(1,x_size); // last grad(E)
  *direction = vector(1,x_size); // descent direction to be taken
  *hessian = vector(1,x_size*x_size); // flattened hessian matrix

  // these vectors are just dummy vectors to use in calculating
  // derivatives, they are allocated once to save time.
  *E_p = vector(1,x_size);
  *E_m = vector(1,x_size);
  *E_pij = vector(1,x_size);
  *E_mij = vector(1,x_size);


  return;
}

void free_vectors(double x_size,double **lastx,double **dEdx,double **lastdEdx,
		  double **direction,double **hessian,double **E_p,
		  double **E_m,double **E_pij,double **E_mij)

{
  free_vector(*lastx,1,x_size);
  free_vector(*dEdx,1,x_size);
  free_vector(*lastdEdx,1,x_size);
  free_vector(*direction,1,x_size);
  free_vector(*hessian,1,x_size*x_size);

  free_vector(*E_p,1,x_size);
  free_vector(*E_m,1,x_size);
  free_vector(*E_pij,1,x_size);
  free_vector(*E_mij,1,x_size);
  return;
}

void allocate_matrices(struct arr_ns ns,double ****c,double ***s,double **r,
		       double ***y,double **r_cp,double ***y_cp,
		       double **rf_fib,int mpt)
/*==============================================================================
  ============================================================================*/
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


void free_matrices(struct arr_ns ns,double ****c,double ***s,double **r,
		   double ***y,double **r_cp,double ***y_cp,double **rf_fib,
		   int mpt)
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

void array_constant(double constant,double *a,int length)
{
  int i;
  for (i = 1; i <= length; i++) a[i] = constant;
  return;
}



void copy_2_arrays(double *r_cp,double **y_cp,double *r,double **y,
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

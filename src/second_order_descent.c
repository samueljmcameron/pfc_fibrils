/* Second order methods of finding the minimum of E(x).

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



bool positive_definite(double *hessian,int x_size)
/*==============================================================================

  Purpose: Determine whether the hessian matrix is positive definite or not by
  computing its eigenvalues.

  ------------------------------------------------------------------------------

  Parameters:

  hessian[1..x_size*x_size] -- the flattened Hessian matrix of E(x).

  x_size -- the number of relevant x vectors in the calculation.

  ------------------------------------------------------------------------------

  Returns: true if the Hessian matrix is positive definite, false otherwise.

  ============================================================================*/
{

  bool positive_eigen(gsl_vector *eigenvals, int x_size);

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,x_size,x_size);
  gsl_matrix *mcp = gsl_matrix_alloc(x_size,x_size);
  gsl_vector *eval = gsl_vector_alloc(x_size);
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(x_size);


  gsl_matrix_memcpy(mcp,&m.matrix);
  gsl_eigen_symm(mcp,eval,w);
  
  int i;

  for (i = 0; i < x_size; i++) {
    printf("eigenvalue_%d = %.6e\n",i,gsl_vector_get(eval,i));
  }
  
  if (!positive_eigen(eval,x_size)) {
    
    gsl_vector_free(eval);
    gsl_eigen_symm_free(w);
    gsl_matrix_free(mcp);
    return false;
  }

  gsl_vector_free(eval);
  gsl_eigen_symm_free(w);
  gsl_matrix_free(mcp);

  return true;
}

bool positive_eigen(gsl_vector *eigenvals, int x_size)
/*==============================================================================

  Purpose: Given a vector of eigenvalues, determine whether they are all
  greater than zero or not.

  ------------------------------------------------------------------------------

  Parameters:

  eigenvals[1..x_size] -- A vector of eigenvalues.

  x_size -- The number of eigenvalues.

  ------------------------------------------------------------------------------

  Returns: true if all the eigenvalues are positive, and false otherwise.

  ============================================================================*/
{
  int i;
  for (i = 0; i < x_size; i++) {
    if (gsl_vector_get(eigenvals,i) <= 0) return false;
  }
  return true;
}

void hessian_update_x(double *x,double *hessian, double *dEdx,int x_size)
/*==============================================================================

  Purpose: Update the vector x[1..x_size] by solving the matrix equation Adx=b,
  where A is the hessian, dx is the change that will be added to x, and b is
  the gradient vector.

  ------------------------------------------------------------------------------

  Parameters:

  x[1..x_size] -- The vector x = (R,eta,delta)' which will be updated.

  hessian[1..x_size*x_size] -- The Hessian of E(x).

  dEdx[1..x_size] -- The gradient of E(x).

  x_size -- The number of relevant quantities in the vector x.

  ------------------------------------------------------------------------------

  Returns: Does not explicitly return anything, but updates the vector x.

  ============================================================================*/
{
  int s;
  int i;

  gsl_vector *dx = gsl_vector_alloc(x_size);

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,x_size,x_size);

  gsl_matrix *mcp = gsl_matrix_alloc(x_size,x_size);

  gsl_vector_view b = gsl_vector_view_array(dEdx+1,x_size);
  gsl_permutation *perm = gsl_permutation_alloc(x_size);

  gsl_matrix_memcpy(mcp,&m.matrix);

  gsl_linalg_LU_decomp(mcp,perm,&s);

  gsl_linalg_LU_solve(mcp,perm,&b.vector,dx);

  for (i = 1; i <= x_size; i++) {
    x[i] -= gsl_vector_get(dx,i-1);
  }

  gsl_permutation_free(perm);
  gsl_matrix_free(mcp);


  return;
}


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


// conjugate gradient descent tools. //

bool positive_eigen(gsl_vector *eigenvals, int x_size);

bool positive_definite(double *hessian, int x_size);

void hessian_update_p(struct params *p, double *hessian, double *dEdx,
		      double *dx,int x_size);


bool positive_definite(double *hessian,int x_size)
{

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
{
  int i;
  for (i = 0; i < x_size; i++) {
    if (gsl_vector_get(eigenvals,i) < 0) return false;
  }
  return true;
}

void hessian_update_p(struct params *p, double *hessian, double *dEdx,
		      double *dx,int x_size)
{

  gsl_vector_view dxcp = gsl_vector_view_array(dx+1,x_size);

  gsl_matrix_view m = gsl_matrix_view_array(hessian+1,x_size,x_size);

  gsl_matrix *mcp = gsl_matrix_alloc(x_size,x_size);
  
  gsl_matrix_memcpy(mcp,&m.matrix);

  gsl_vector_view b = gsl_vector_view_array(dEdx+1,x_size);
  int s;
  gsl_permutation *perm = gsl_permutation_alloc(x_size);

  //  compute_eigenvalues(hessian);

  gsl_linalg_LU_decomp(mcp,perm,&s);

  gsl_linalg_LU_solve(mcp,perm,&b.vector,&dxcp.vector);

  p->R -= gsl_vector_get(&dxcp.vector,0);
  if (p->delta > EFFECTIVE_ZERO) {
    p->eta -= gsl_vector_get(&dxcp.vector,1);
    p->delta -= gsl_vector_get(&dxcp.vector,2);
  }

  gsl_permutation_free(perm);
  gsl_matrix_free(mcp);


  return;
}


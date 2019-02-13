/* finite difference calculations for first order and second order derivatives */


#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>
#include "nrutil.h"
#include "headerfile.h"




void c_deriv(double (*f)(const gsl_vector *,void *),gsl_vector *x,int i,
	     void *ps,double h,double *result,double *abserr_round,
	     double *abserr_trunc)
{
  /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
     x+h/2, x+h). Note that the central point is not used.  

     Compute the error using the difference between the 5-point and
     the 3-point rule (x-h,x,x+h). Again the central point is not
     used. */

  gsl_vector_set(x,i,gsl_vector_get(x,i)-h);
  double fm1 = f(x,ps);
  gsl_vector_set(x,i,gsl_vector_get(x,i)+h+h);
  double fp1 = f(x,ps);

  gsl_vector_set(x,i,gsl_vector_get(x,i)-h-h/2.0);
  double fmh = f(x,ps);
  gsl_vector_set(x,i,gsl_vector_get(x,i)+h);
  double fph = f(x,ps);
  gsl_vector_set(x,i,gsl_vector_get(x,i)-h/2.0);

  double r3 = 0.5 * (fp1 - fm1);
  double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

  double e3 = (fabs (fp1) + fabs (fm1)) * GSL_DBL_EPSILON;
  double e5 = 2.0 * (fabs (fph) + fabs (fmh)) * GSL_DBL_EPSILON + e3;

  /* The next term is due to finite precision in x+h = O (eps * x) */

  double dy = (GSL_MAX(fabs(r3/h),fabs(r5/h))
	       *(fabs(gsl_vector_get(x,i))/ h)*GSL_DBL_EPSILON);

  /* The truncation error in the r5 approximation itself is O(h^4).
     However, for safety, we estimate the error from r5-r3, which is
     O(h^2).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r5. */

  *result = r5 / h;
  *abserr_trunc = fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
  *abserr_round = fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
  return;
}

int deriv_xi(double (*f)(const gsl_vector *,void *),const gsl_vector *x_scale,
	     int i,void *ps,double h,double *result,double *abserr)
{

  void c_deriv(double (*f)(const gsl_vector *,void *),gsl_vector *x,int i,
	       void *ps,double h,double *result,double *abserr_round,
	       double *abserr_trunc);
  
  
  double r_0, round, trunc, error;
  gsl_vector *x;

  struct params *p = ps;

  x=gsl_vector_alloc(p->x_size);
  gsl_vector_memcpy(x,x_scale);
  c_deriv(f,x,i,ps,h,&r_0,&round,&trunc);
  error = round + trunc;

  if (round < trunc && (round > 0 && trunc > 0))
    {
      double r_opt, round_opt, trunc_opt, error_opt;

      /* Compute an optimised stepsize to minimize the total error,
         using the scaling of the truncation error (O(h^2)) and
         rounding error (O(1/h)). */

      double h_opt = h * pow (round / (2.0 * trunc), 1.0 / 3.0);
      c_deriv(f,x,i,ps,h_opt,&r_opt,&round_opt,&trunc_opt);
      error_opt = round_opt + trunc_opt;

      /* Check that the new error is smaller, and that the new derivative 
         is consistent with the error bounds of the original estimate. */

      if (error_opt < error && fabs (r_opt - r_0) < 4.0 * error) {
	r_0 = r_opt;
	error = error_opt;
      }
    }

  *result = r_0;
  *abserr = error;

  return GSL_SUCCESS;

}

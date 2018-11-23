#ifndef _HEADERFILE_H_
#define _HEADERFILE_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_vector.h>

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#define NE 2                    // # of 1st order DEs
#define MAX_M ((256)*(256)+(1)) // max # of mesh points (2^M+1 for romberg integration)
#define NB 1                    // # of BCs at first boundary (k = 1)
#define NSI NE                  // max # i of S_i,j
#define NSJ (2*NE+1)            // max # j of S_i,j
#define NYJ NE                  // # of dependent variables
#define NYK M                   // number of points in each dependent variable
#define NCI NE                  // # of rows in storage matrix within c[m][:][:]
#define NCJ (NE-NB+1)           // # of columns in storage matrix within c[m][:][:]
#define NCK (M+1)               // # number of points in tensor c, c[:][m][n]

#define X_SIZE     3
#define ITMAX      5000
#define CONV_ODE   1e-10
#define CONV_MIN   1e-7



struct params{
  double K33;
  double k24;
  double Lambda;
  double d0;
  double omega;
  double gamma_s;
  double Rguess;
  double etaguess;
  double deltaguess;
  double Rupper;
  double Rlower;
  double etaupper;
  double etalower;
  double deltaupper;
  double deltalower;
  double upperbound_x;
  double upperbound_y;
  double *r;
  double **y;
  double *rf_fib;
  double *r_cp;
  double **y_cp;
  double **s;
  double ***c;
  int mpt;
};


/* from file energy.c */

double E_calcwrap(double *x_scale,void *ps);

/* from file scaling.c */

void scale_forward(gsl_vector *y,const double *x,struct params *p);
void scale_backward(const gsl_vector *y, double *x,struct params *p);


#endif /* ANSI */

#endif /* _HEADERFILE_H_ */
#ifndef _PROJECTILE_H_
#define _PROJECTILE_H_

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

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */


struct arr_ns{
  int ne; // number of equations in ODE (2)
  int nb; // number of BCs at first boundary
  int nsi;
  int nsj;
  int nyj;
  int nci;
  int ncj;
};

struct params{
  double K33;
  double k24;
  double Lambda;
  double d0;
  double omega;
  double R;
  double eta;
  double delta;
  double gamma_s;
  double upperbound_x;
  double upperbound_y;
};

void solvde(int itmax, double conv, double slowc, double scalv[],
	    struct arr_ns *ns, int m, double **y, double *r,
	    double ***c, double **s, struct params *p,
	    double h);

void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,
	   struct params *p, double h, int mpt);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s);
double trapzd(double *, double *, double, int, int);
void polint(double xa[], double ya[], int, double, double *, double *);
double qromb(double *,double *, int,double tol,bool *failure);

bool energy_stuff(double *E, double *dEdR,double *dEdeta,
		  double *dEddelta,struct params *p,
		  double *r,double **y,double *rf_,
		  double *integrand1,double *integrand2,int mpt);

void scanE(struct params p,FILE *energy,FILE *psi,double conv,
	   int itmax,int mpt,int num_x,char scan_what[]);

void scan2dE(struct params p,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,int num_x, int num_y,
	     char scan_what_x[],char scan_what_y[]);

void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,
	      FILE *denergyddelta,FILE *surfacetwist,
	      double conv,int itmax,int mpt,double rate);


#endif /* ANSI */

#endif /* _PROJECTILE_H_ */

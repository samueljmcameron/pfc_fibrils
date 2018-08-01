#ifndef _PROJECTILE_H_
#define _PROJECTILE_H_

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

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void solvde(int itmax, double conv, double slowc, double scalv[],
	    int ne, int nb, int m, double **y, double *r,
	    double ***c, double **s, double K33, double k24,
	    double Lambda, double eta, double d0, double L,
	    double h,int mpt);
void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,double K33,
	   double k24, double Lambda, double eta, double d0,
	   double L,double h, int mpt);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s);
double trapzd(double *, double *, double, int, int);
void polint(double xa[], double ya[], int, double, double *, double *);
double qromb(double *,double *, int);

void compute_rf2233b(double K33, double Lambda,double eta, 
		     double d0, double L,double *r,
		     double **y, double *rf_,int mpt);

void compute_integrand1(double Lambda,double eta, double d0,
			double L, double *r, double **y,
			double *integrand1,int mpt);

void compute_integrand2(double Lambda,double eta, double d0,
			double L, double *r, double **y,
			double *integrand2,int mpt);

double E_R(double R, double k24, double L, double gamma_s,
	   double gamma_t,double *r,double **y, 
	   double integration_2233b,int mpt);

double derivEdR(double R, double K33, double k24,double Lambda,
		double eta, double d0, double L, double gamma_s,
		double *r,double **y, double integration_2233b,
		int mpt);

double derivEdeta(double R,double Lambda,double eta, 
	      double d0, double L, double *r,double **y,
	      double integration1,double integration2);

double derivEdL(double R,double Lambda,double eta, 
	    double d0, double L, double gamma_t,
	    double *r,double **y,double integration2);

void energy_stuff(double *E, double *dEdR,double *dEdeta, double *dEdL,
		  double R, double K33, double k24, double Lambda,
		  double eta, double d0, double L, double gamma_s,
		  double gamma_t, double *r, double **y, double *rf_,
		  double *integrand1, double *integrand2,int mpt);


#endif /* ANSI */

#endif /* _PROJECTILE_H_ */

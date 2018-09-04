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

void solvde(int itmax, double conv, double slowc, double scalv[],
	    int ne, int nb, int m, double **y, double *r,
	    double ***c, double **s, double K33, double k24,
	    double Lambda,double d0,double eta,double delta,
	    double h);
void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,double K33,
	   double k24, double Lambda, double d0,double eta,
	   double delta, double h, int mpt);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s);
double trapzd(double *, double *, double, int, int);
void polint(double xa[], double ya[], int, double, double *, double *);
double qromb(double *,double *, int,bool *failure);

void compute_rf2233b1(double K33, double Lambda,double d0,
		      double eta,double delta,double *r,
		      double **y, double *rf_,int mpt);
void compute_integrand1(double d0,double eta,double *r,
			double **y,double *integrand1,int mpt);
void compute_integrand2(double d0,double eta,double *r,
			double **y,double *integrand2,int mpt);


double E_R(double k24,double Lambda,double omega,double R,
	   double eta,double delta,double gamma_s,double *r,
	   double **y,double integration_2233b1,int mpt);

double derivEdR(double K33, double k24,double Lambda,double d0, 
		double R,double eta,double delta,double gamma_s,
		double *r,double **y,double integration_2233b1,
		int mpt);

//double derivEdL(double Lambda,double omega,double R,double eta,
//		double delta, double gamma_t,double integration2);

double derivEdeta(double Lambda,double omega,double R,double eta,
		  double delta,double integration1,double integration2);

double derivEddelta(double Lambda,double omega,double R,double eta,
		    double delta,double integration2);

bool energy_stuff(double *E, double *dEdR,double *dEdeta,
		  double *dEddelta,double K33, double k24,
		  double Lambda,double d0,double omega,double R,
		  double eta,double delta,double gamma_s,
		  double *r,double **y,double *rf_,
		  double *integrand1,double *integrand2,int mpt);

void write_failure(double *r, double **y,double *rf_,int rlength,char *f_err);

void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt);

void propagate_r(double *r, double h,int mpt);

void save_psi(FILE *psi,double *r, double **y,int mpt);

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable);

void make_f_err(char *f_err,int f_err_size,double K33,double k24,
		double Lambda,double d0,double omega,double R,
		double eta,double delta,double gamma_s);

void scanE(double K33,double k24,double Lambda,double d0,
	   double omega,double R,double eta,double delta,
	   double gamma_s,FILE *energy,FILE *psi,
	   double conv,int itmax,int mpt,
	   double upperbound,char scan_what[]);

void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_xpoints);

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int xpoints);

void quick_interp(double *xa, double *ya, double x, double *y,
		  int xpoints);

void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],double *R, 
			double *dEdR,double *dEdRlast,double *eta,
			double *dEdeta, double *dEdetalast,double *delta,
			double *dEddelta,double *dEddeltalast);

void scan2dE(double *r,double **y,double ***c,double **s,
	     double K33,double k24,double Lambda,double d0,
	     double omega,double R,double eta,double delta,
	     double gamma_s,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,double upperbound_x,double upperbound_y,
	     char scan_what_x[],char scan_what_y[]);

void resize_all_arrays(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int xpoints,int nci,int ncj,
		       int nsi,int nsj,int nyj);

#endif /* ANSI */

#endif /* _PROJECTILE_H_ */

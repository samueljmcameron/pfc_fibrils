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
  //  double R;
  //  double eta;
  //  double delta;
  double gamma_s;
  double upperbound_x;
  double upperbound_y;
};

bool solvde(int itmax, double conv, double slowc, double scalv[],
	    struct arr_ns *ns, int m, double **y, double *r,
	    double ***c, double **s, struct params *p,
	    double h);

void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,
	   struct params *p, double h, int mpt);
bool pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s);
double trapzd(double *, double *, double, int, int);
void polint(double xa[], double ya[], int, double, double *, double *);
double qromb(double *,double *, int,double tol,bool *failure);


bool energy_properties(double *E,double *dEdx,struct params *p,
		       double *r,double **y,double *rf_,
		       double *integrand1,double *integrand2,int mpt);

bool energy_prop_with_hessian(double *E, double *dEdx,struct params *p,
			      double *r,double **y,double *rf_,
			      double *integrand1,double *integrand2,int mpt,
			      double *hessian);

void scanE(struct params p,FILE *energy,FILE *psi,double conv,
	   int itmax,int mpt,int num_scan,char scan_what[]);

void scan2dE(struct params p,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,int num_scanx, int num_scany,
	     char scan_what_x[],char scan_what_y[]);

void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,
	      FILE *denergyddelta,FILE *surfacetwist,
	      FILE *energydensity,double conv,int itmax,
	      int mpt,double rate0);


/* files below this point are in shared.c */

// energy (optional hessian) calculation function (links this file to energy.c file). //

void single_calc(double *E,double *dEdx,struct params *p,
		 double ****c,double ***s,double ***y,double **r,
		 double **rf_,double **integrand1,double **integrand2,
		 double **y_cp,double *r_cp,double *hessian,
		 double conv,int itmax,int *npoints,int last_npoints,
		 struct arr_ns *ns,int max_size,bool calc_Hess);

// file I/O functions.

void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p);

void write_failure(char *err_type,double *r, double **y,double *rf_,
		   int rlength,struct params p);

void save_psi(FILE *psi,double *r, double **y,int mpt);

void save_energydensity(FILE *energydensity,double *r, double *rf_,int mpt);

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable);



// array and utility functions which are specific to algorithms being used. //

void assign_ns(struct arr_ns *ns);

void resize_and_interp(double h,double ****c,double ***s,double ***y,double **r,
		       double **rf_,double **integrand1,double **integrand2,
		       double **y_cp,double *r_cp,int npoints,int last_npoints,
		       struct arr_ns *ns);

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns);

void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns);

void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt);

void propagate_r(double *r, double h,int mpt);

void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],
			struct params *p,double *dEdx,double *lastdEdx);



// array interpolation utility functions. //

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int npoints);

void quick_interp(double *xcp,double **ycp,double x,double **y,int i);


// utility functions that are not specific to algorithms being used. //

int sign(double x);

void arr_cp(double *acp, double *a,int length);

double vector_norm(double *a,int length);

void array_constant(double constant,double *a,int length);

bool non_zero_array(double *dEdx,double conv);

void copy_2_arrays(double *r,double **y,double *r_cp,double **y_cp,
		   int last_npoints);


#endif /* ANSI */

#endif /* _PROJECTILE_H_ */

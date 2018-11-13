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

/* files from shared.c */

void assign_ns(struct arr_ns *ns);

void allocate_matrices(struct arr_ns ns,double ****c,double ***s,double **r,
		       double ***y,double **r_cp,double ***y_cp,
		       double **rf_fib,int mpt);

void allocate_vectors(double x_size,double **lastx,double **dEdx,
		      double **lastdEdx,double **direction,double **hessian,
		      double **E_p,double **E_m,double **E_pij,double **E_mij);

void free_vectors(double x_size,double **lastx,double **dEdx,double **lastdEdx,
		  double **direction,double **hessian,double **E_p,
		  double **E_m,double **E_pij,double **E_mij);

void free_matrices(struct arr_ns ns,double ****c,double ***s,double **r,
		   double ***y,double **r_cp,double ***y_cp,double **rf_fib,
		   int mpt);

void array_constant(double constant,double *a,int length);

void arr_cp(double *acp, double *a,int length);

void copy_2_arrays(double *r,double **y,double *r_cp,double **y_cp,
		   int last_mpt);

void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p,
		double *x);

void save_psi(FILE *psi,double *r, double **y,int mpt);

void save_energydensity(FILE *energydensity,double *r, double *rf_fib,
			int mpt);

bool non_zero_array(double *dEdx,double conv,int x_size);

/* from file solvde.c */

void linearGuess(double *r, double **y, double initialSlope,double h,int mpt);

void solvde_wrapper(int itmax, double convODE, double scalv[],struct arr_ns *ns,
		    int mpt,double *r,double **y,double **y_guess,double ***c,
		    double **s,struct params *p,double *x,double h);

bool solvde(int itmax, double conv, double scalv[],struct arr_ns *ns, int m,
	    double *r, double **y,double ***c, double **s,struct params *p,
	    double *x,double h);

/* from file energy.c */

double E_calc(struct params *p,double *x,double *r,double **y,double *rf_fib,
	      double ***c,double **s,double *r_cp,double **y_cp,double convODE,
	      int itmax,int *mpt,struct arr_ns *ns,int max_mpt);

bool successful_E_count(double *E,struct params *p,double *x,double *r,
			double **y,double *rf_fib,int mpt);

/* from file finite_differences.c */

void derivatives_fd(double *dEdx,double E,struct params *p,double *x,
		    double ***c,double **s,double *r,double **y,double *rf_fib,
		    double *r_cp,double **y_cp,double convODE,double dx,
		    int itmax,int *mpt,struct arr_ns *ns,int max_mpt,
		    int x_size,double *hessian,bool calc_hess,double *E_p,
		    double *E_m, double *E_pij,double *E_mij);

double compute_dx(double E,double convMIN);


/* from file conj_grad.c */

void set_direction(double *direction,double *dEdx,double *lastdEdx,int x_size);

void armijo_backtracker(double rate,double E,double *dEdx,double *direction,
			struct params *p,double *x,double *r,double **y,
			double *rf_fib, double ***c, double **s,double *r_cp,
			double **y_cp,double convODE,int itmax,int *mpt,
			struct arr_ns *ns,int max_mpt,double min_rate,
			int x_size);

/* from second_order_descent.c */

bool positive_definite(double *hessian,int x_size);

void hessian_update_x(double *x,double *hessian, double *dEdx,int x_size);

/* from graddesc_driver.c */

void graddesc(struct params p,double *x,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,FILE *denergyddelta,
	      FILE *surfacetwist,FILE *energydensity,const double convODE,
	      const double convMIN,const int itmax,int mpt,const int max_mpt,
	      double rate,const int x_size0);

/* from shooting.c */
void shootsolve_driver(struct params p, double *x, FILE *psivsr,double psip01,
		       double psip02,int mpt);

void shootscan_driver(struct params p, double *x, FILE *bc,double psip01,
		      double psip02,int numpoints,int mpt);


/* functions from Numerical Recipes, in files matching function names. */

void bksub(int ne, int nb, int jf, int k1, int k2, double ***c);
void difeq(int k, int k1, int k2, int jsf, int isl, int isf,
	   int ne, double **s, double **y, double *r,
	   struct params *p, double *x,double h, int mpt);
bool pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s);
double trapzd(double *, double *, double, int, int);
void polint(double xa[], double ya[], int, double, double *, double *);
double qromb(double *,double *, int,double tol,bool *failure);




#endif /* ANSI */

#endif /* _PROJECTILE_H_ */

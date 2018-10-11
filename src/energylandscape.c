#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"
#include <time.h>

int sign(double x) {
  return (x > 0) - (x < 0);
}

void arr_cp(double *acp, double *a,int length)
{
  int i;
  for (i = 1; i <= length; i++) acp[i] = a[i];
  return;
}

double vector_norm(double *a,int length)
{
  int i;
  double sum = 0;
  for (i = 1; i <= length; i++) sum += a[i]*a[i];

  return sum;
}

void array_constant(double constant,double *a,int length)
{
  int i;
  for (i = 1; i <= length; i++) a[i] = constant;
  return;
}

bool non_zero_array(double *dEdx,double conv)
{
  return (fabs(dEdx[1]) > conv || fabs(dEdx[2]) > conv
	  || fabs(dEdx[3]) > conv);
}

void update_p(struct params *p,double rate,double *arrchange)
{
  p->R -= rate*arrchange[1];
  p->eta -= rate*arrchange[2];
  p->delta -= rate*arrchange[3];
  return;
}

void reset_p(struct params *p,double rate,double *arrchange)
{
  p->R += rate*arrchange[1];
  p->eta += rate*arrchange[2];
  p->delta += rate*arrchange[3];
  return;
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
void write_failure(char *err_type,double *r, double **y,double *rf_,
		   int rlength,struct params p);


void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt);
void propagate_r(double *r, double h,int mpt);
void save_psi(FILE *psi,double *r, double **y,int mpt);
void save_energydensity(FILE *energydensity,double *r, double *rf_,int mpt);
void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable);
void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p);
void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_npoints);
void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int npoints);
void quick_interp(double *xcp,double **ycp,double x,double **y,int i);
void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],
			struct params *p,double *dEdx,double *lastdEdx);

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns);
void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns);
void single_calc(double *E,double *dEdx,struct params *p,
		 double ****c,double ***s,double ***y,double **r,
		 double **rf_,double **integrand1,double **integrand2,
		 double **y_cp,double *r_cp,double conv,int itmax,
		 int *npoints,int last_npoints,struct arr_ns *ns,
		 int max_size);
double backtracker(double beta,double rate,double E,double *dEdx,
		   struct params *p,double *r,double **y,double *r_cp,
		   double **y_cp,double conv, int itmax,int npoints,
		   struct arr_ns *ns,int last_npoints,int max_size,
		   double min_rate);
bool jumpmin(double frac_tol,double E,double *dEdR,struct params p,
	     double ***c,double **s,double **y,double *r,double *rf_,
	     double *integrand1,double *integrand2,double **y_cp,
	     double *r_cp,double conv,int itmax,int npoints,
	     int last_npoints,struct arr_ns *ns,int max_size);

void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt)
{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess!
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*r[k]; // y1 is psi
    y[2][k] = initialSlope; // y2 is psi'!!!!!!!!!!!!!!!!!!
  }
  return;
}

void propagate_r(double *r, double h,int mpt)
{
  int k;
  for (k=1;k <=mpt; k++) r[k] = (k-1)*h; // only change r since psi, psi' are stored from last loop
  return;
}

void save_psi(FILE *psi,double *r, double **y,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(psi,"%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i]);
  }
  printf("psi(R) = %1.2e\n",y[1][mpt]);
  return;
}

void save_energydensity(FILE *energydensity,double *r, double *rf_,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(energydensity,"%10.8e\t%10.8e\n",r[i],rf_[i]);
  }
  return;
}

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable)
{
  fprintf(energy,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,E,
	  derivative,observable);
  return;
}

void make_f_err(char *f_err,char *err_type,int f_err_size,struct params p)
{
  snprintf(f_err,f_err_size,"data/%s_psivsr_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",err_type,
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,
	   p.delta,p.gamma_s);
  return;
}


void copy_arrays(double *r,double **y,double *r_cp,double **y_cp,
		 int last_npoints)
{
  int i;

  for (i = 1; i <=last_npoints; i++) {
    y_cp[1][i] = y[1][i];
    y_cp[2][i] = y[2][i];
    r_cp[i] = r[i];
  }

  return;
}

void interpolate_array(double *r,double **y,double *r_cp,
		       double **y_cp,int npoints)
{
  int i;
  double dy;

  for (i = 1; i <=npoints; i++) quick_interp(r_cp,y_cp,r[i],y,i);
  return;
}

void quick_interp(double *xcp,double **ycp,double x,double **y,int i)
// Given two points of ycp[1][:] (y11cp,x1cp) and y(12cp,x2cp), //
// and two points of ycp[2][:] (y21cp,x1cp) and (y22cp,x2cp),   //
// interpolate to determine the value of y between the two

{
  if (i%2 != 0) {
    y[1][i] = ycp[1][(i-1)/2+1];
    y[2][i] = ycp[2][(i-1)/2+1];
  }
  else {
    double slope_1, slope_2;
    slope_1 = (ycp[1][i/2+1]-ycp[1][i/2])/(xcp[i/2+1]-xcp[i/2]);
    slope_2 = (ycp[2][i/2+1]-ycp[2][i/2])/(xcp[i/2+1]-xcp[i/2]);

    y[1][i] = slope_1*(x-xcp[i/2])+ycp[1][i/2];
    y[2][i] = slope_2*(x-xcp[i/2])+ycp[2][i/2];
  }

  return;
}

void setup_var_pointers(double **var, double *var0,double **dEdvar,
			double **dEdvarlast,char scan_what[],
			struct params *p,double *dEdx,double *lastdEdx)
// convenient way to make a variables (ending in var) that have the  //
// same addresses as whichever parameter that is specified by        //
// scan_what.
{

  if (strcmp(scan_what,"R")==0) {
    printf("R!\n");
    *var = &p->R;
    *var0 = p->R;
    *dEdvar = &dEdx[1];
    *dEdvarlast = &lastdEdx[1];
  }
  else if (strcmp(scan_what,"eta")==0) {
    printf("eta!\n");
    *var = &p->eta;
    *var0 = p->eta;
    *dEdvar = &dEdx[2];
    *dEdvarlast = &lastdEdx[2];
  }
  else if (strcmp(scan_what,"delta")==0) {
    printf("delta!\n");
    *var = &p->delta;
    *var0 = p->delta;
    *dEdvar = &dEdx[3];
    *dEdvarlast = &lastdEdx[3];
  }
  else {
    printf("Need either R, eta, or delta as argv[n] input."
	   "Exiting to system.\n");
    exit(1);
  }
  return;
}


void write_failure(char *err_type,double *r, double **y,double *rf_,int rlength,struct params p)
{
  int i;
  FILE *broken;
  int f_err_size = 200;
  char f_err[f_err_size];

  make_f_err(f_err,err_type,f_err_size,p);

  if (strcmp(err_type,"QROMB")==0) {
      printf("failed to integrate at npoints = %d, with R = %e.\n",rlength,r[rlength]);
      printf("saving psi(r) shape, and exiting to system.\n");
    }
  broken = fopen(f_err,"w");


  for (i = 1; i<=rlength; i++) {
    fprintf(broken,"%.8e\t%.8e\t%.8e\t%.8e\n",r[i],y[1][i],y[2][i],rf_[i]);
  }
  fclose(broken);
  exit(1);
  return;
}

void resize_all_arrays(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns)
{
  int last_npoints;

  last_npoints = (npoints-1)/2+1;
  free_matrices(c,s,y,r,rf_,integrand1,integrand2,last_npoints,ns);
  allocate_matrices(c,s,y,r,rf_,integrand1,integrand2,npoints,ns);


  return;
}

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns)
{
  *y = matrix(1,ns->nyj,1,npoints);
  *s = matrix(1,ns->nsi,1,ns->nsj);
  *c = f3tensor(1,ns->nci,1,ns->ncj,1,npoints+1);
  *r = vector(1,npoints);
  *rf_ = vector(1,npoints);
  *integrand1 = vector(1,npoints);
  *integrand2 = vector(1,npoints);
  return;
}

void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns)
{

  free_f3tensor(*c,1,ns->nci,1,ns->ncj,1,npoints+1);
  free_matrix(*s,1,ns->nsi,1,ns->nsj);
  free_matrix(*y,1,ns->nyj,1,npoints);
  free_vector(*r,1,npoints);
  free_vector(*rf_,1,npoints);
  free_vector(*integrand1,1,npoints);
  free_vector(*integrand2,1,npoints);
  return;
}


void scanE(struct params p,FILE *energy,FILE *psi,double conv,
	   int itmax,int mpt,int num_x,char scan_what[])
// The energy of the system E(R,L,eta). This function  //
// generates data of E vs scan_what[] (either "delta", //
// "R", or "eta"), while holding the other two values  //
// constant. The E vs scan_what[] data is saved in     //
// energy file. If a minimum (or multiple minima) are  //
// found in the energy landscape, psi(r) is saved to   //
// the psi file. //
{
  int npoints = mpt;
  int last_npoints = mpt;
  double *var,var0;
  double h;
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double *dEdx, *lastdEdx;
  double initialSlope;
  double E;
  double *dEdvar, *dEdvarlast;
  double Emin = 1e100;
  int max_size = (mpt-1)*4+1;
  double dx;
  int count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);

  setup_var_pointers(&var,&var0,&dEdvar,&dEdvarlast,scan_what,
		     &p,dEdx,lastdEdx);

  dx = (p.upperbound_x-var0)/num_x;
  printf("dx = %lf\n",dx);
  
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  initialSlope = M_PI/(4.0*p.R);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  count_x = 0;

  while (count_x <= num_x) {

    // for each value of x (i.e *var) in E vs x


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);

    last_npoints = npoints;

    // save var,E, and surface twist
    saveEnergy(energy,*var,E,*dEdvar,y[1][npoints]);
      


    // if a minimum has been found, save the psi(r) curve
    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0
	&& *dEdvarlast < 0 && E <= Emin) {
      save_psi(psi,r,y,npoints);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E);
      Emin = E;
    }

    count_x += 1;
    *var = var0+dx*count_x;
    arr_cp(lastdEdx,dEdx,3);
    
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);

  return;
}



void scan2dE(struct params p,FILE *energy,FILE *psi,
	     FILE *deriv_energy_x,FILE *deriv_energy_y,
	     FILE *surfacetwist,double conv,int itmax,
	     int mpt,int num_x, int num_y,
	     char scan_what_x[],char scan_what_y[])
// The energy of the system E(R,L,eta). This function    //
// generates data of E vs scan_what_x[] vs scan_what_y[] //
// (either "delta","R", or "eta"), while holding the     //
// other value constant. The E vs x vs y data is saved   //
// in energy file. If a minimum (or multiple minima) are //
// found in the energy landscape, psi(r) are saved to    //
// the psi file. //
{
  int npoints = mpt;
  int last_npoints = mpt;
  double *var_x,var_x0;
  double *var_y,var_y0;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double *dEdx, *lastdEdx;
  double *dEdvar_x, *dEdvar_xlast;
  double *dEdvar_y, *dEdvar_ylast;
  double Emin = 1e100;
  int max_size = (mpt-1)*10+1;
  double dx,dy;
  int count_y,count_x;
  struct arr_ns ns;

  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);

  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  

  setup_var_pointers(&var_x,&var_x0,&dEdvar_x,&dEdvar_xlast,scan_what_x,
		     &p,dEdx,lastdEdx);
  dx = (p.upperbound_x-var_x0)/(num_x-1);

  setup_var_pointers(&var_y,&var_y0,&dEdvar_y,&dEdvar_ylast,scan_what_y,
		     &p,dEdx,lastdEdx);
  dy = (p.upperbound_y-var_y0)/(num_y-1);

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  printf("dx = %lf, dy = %lf\n",dx,dy);
  
  initialSlope = M_PI/(4.0*p.R);

  count_y = 0;

  while (count_y < num_y) {
    npoints = mpt;
    last_npoints = mpt;
    count_x = 0;
    *var_x = var_x0;
    h = p.R/(npoints-1);
    // initial guess for the functional form of psi(r) and psi'(r)
    printf("initial slope for guess = %lf\n", initialSlope);
    linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

    
    while (count_x < num_x) {

      single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		  y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);

      last_npoints = npoints;

      if (count_x == 0) initialSlope = 0.1*y[2][npoints];

      fprintf(energy,"%.8e\t",E);
      fprintf(deriv_energy_x,"%.8e\t",*dEdvar_x);
      fprintf(deriv_energy_y,"%.8e\t",*dEdvar_y);
      fprintf(surfacetwist,"%.8e\t",y[1][mpt]);

      if (*var_x != var_x0 && *dEdvar_x*(*dEdvar_xlast) <= 0
	  && *dEdvar_xlast < 0 && *var_y != var_y0
	  && *dEdvar_y*(*dEdvar_ylast) <= 0 && *dEdvar_ylast <0) {
	save_psi(psi,r,y,npoints);
	printf("SAVED!\n");
	printf("E_min-E_chol = %1.2e\n",E);
	Emin = E;
      }

      arr_cp(lastdEdx,dEdx,3);
      count_x += 1;
      *var_x = var_x0+count_x*dx;
    }
    fprintf(energy,"\n");
    fprintf(deriv_energy_x,"\n");
    fprintf(deriv_energy_y,"\n");
    fprintf(surfacetwist,"\n");
    count_y += 1;
    *var_y = var_y0+count_y*dy;
    printf("%s = %lf\n",scan_what_x,*var_x-dx);
    printf("%s = %lf\n",scan_what_y,*var_y-dy);
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);
  
  return;
}

void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,
	      FILE *denergyddelta,FILE *surfacetwist,
	      FILE *energydensity,double conv,int itmax,
	      int mpt,double rate0)
{
  int npoints = mpt;
  int last_npoints = mpt;
  double h;
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double *dEdx, *lastdEdx;
  double beta = 0.7;
  double Elast = -1e100; // a really large, negative number to start with
  double min_rate = 1e3*conv;
  double rate = rate0;
  double lastR = 1e100;
  double lasteta = 1e100;
  double lastdelta = 1e100;
  int max_size = (mpt-1)*8+1;
  int count = 0;
  double frac_tol = 1e-6;
  clock_t begin = clock();
  clock_t end;

  struct arr_ns ns;
  assign_ns(&ns);
  dEdx = vector(1,3);
  lastdEdx = vector(1,3);  
 

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  
  initialSlope = M_PI/(4.0*p.R);
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  printf("initial slope for guess = %lf\n", initialSlope);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  // using classical gradient descent, try to find minimum.


  array_constant(1e100,dEdx,3);
  arr_cp(lastdEdx,dEdx,3);

  while (non_zero_array(dEdx,conv)) {


    single_calc(&E,dEdx,&p,&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		y_cp,r_cp,conv,itmax,&npoints,last_npoints,&ns,max_size);


    last_npoints = npoints;

    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.8e\t%.8e\n",p.R,dEdx[1]);
    fprintf(denergydeta,"%.8e\t%.8e\n",p.eta,dEdx[2]);
    fprintf(denergyddelta,"%.8e\t%.8e\n",p.delta,dEdx[3]);
    fprintf(surfacetwist,"%.8e\t%.8e\n",p.R,y[1][mpt]);

    rate = backtracker(beta,rate0,E,dEdx,&p,r,y,r_cp,y_cp,
		       conv,itmax,npoints,&ns,last_npoints,max_size,min_rate);

    update_p(&p,rate,dEdx);

    count += 1;
    //    printf("count = %d\n",count);

    if (rate <= min_rate) {
      //      printf("last rate <= min_rate, where min_rate = %e and last "
      //     " rate = %e.\n",min_rate,rate);
      if (fabs(lastR-p.R)<conv && fabs(lasteta-p.eta)<conv
	  && fabs(lastdelta-p.delta)<conv && fabs(Elast-E)<conv) {
	//printf("changes in R, eta, and delta are smaller than %e.\n",conv);
	if(jumpmin(frac_tol,E,dEdx,p,c,s,y,r,rf_,integrand1,
		   integrand2,y_cp,r_cp,conv,itmax,npoints,
		   last_npoints,&ns,max_size)) break;
      }
    }

    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    arr_cp(lastdEdx,dEdx,3);
    Elast = E;

    if (count == 10000) {
      end = clock();
      printf("Elapsed after %d: %lf seconds\n",count,
	     (double)(end - begin) / CLOCKS_PER_SEC);
      exit(1);
    }

    if (p.R <= 0) {
      printf("R is being driven to negative values, R = %e!\n",p.R);
      p.R = 1e-6;
    }
      
  }



  



  printf("count = %d\n",count);
  save_psi(psi,r,y,npoints);
  save_energydensity(energydensity,r,rf_,npoints);
  printf("SAVED!\n");
  printf("E_min-E_chol = %1.2e\n",E);

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);
  free_vector(dEdx,1,3);
  free_vector(lastdEdx,1,3);

  
  return;
}

bool jumpmin(double frac_tol,double E,double *dEdx,struct params p,
	     double ***c,double **s,double **y,double *r,double *rf_,
	     double *integrand1,double *integrand2,double **y_cp,
	     double *r_cp,double conv,int itmax,int npoints,
	     int last_npoints,struct arr_ns *ns,int max_size)
{
  int i;
  double lastR = p.R;
  double lasteta = p.eta;
  double lastdelta = p.delta;
  double *tmpdEdx;
  double Elast = E;
  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;

  tmpdEdx = vector(1,3);
  arr_cp(tmpdEdx,dEdx,3);


  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_arrays(r,y,rc,yc,last_npoints);


  p.R = lastR*(1-sign(dEdx[1])*frac_tol);
  p.eta = lasteta*(1-sign(dEdx[2])*frac_tol);
  p.delta = lastdelta*(1-sign(dEdx[3])*frac_tol);

  single_calc(&E,tmpdEdx,&p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,conv,itmax,&npoints,last_npoints,
	      ns,max_size);

  printf("lastE = %e, E = %e\n",Elast,E);
  for (i = 1; i <= 3; i++) {
    printf("initial dEdx[%d] = %e, new dEdx[%d] = %e,",
	   i,dEdx[i],i,tmpdEdx[i]);
  }
  printf("\n\n");
  if (dEdx[1]*tmpdEdx[1] <=0 && dEdx[2]*tmpdEdx[2] <= 0
      && dEdx[3]*tmpdEdx[3] <= 0) {
    printf("jumped over the minimimum (derivatives all changed signs).\n");
    return true;
  }


  return false;
}



void single_calc(double *E,double *dEdx,struct params *p,
		 double ****c,double ***s,double ***y,double **r,
		 double **rf_,double **integrand1,double **integrand2,
		 double **y_cp,double *r_cp,double conv,int itmax,
		 int *npoints,int last_npoints,struct arr_ns *ns,
		 int max_size)
{
  double h;
  bool successful_qromb=false;
  bool successful_solvde;
  int solvde_count;
  double slopeguess = M_PI/(4.0*p->R);
  double slowc = 1.0;
  double scalv[2+1];

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values

  do {
    h = p->R/((*npoints)-1);
    
    // only executes if the first calculation for the current variable values is unsuccessful
    if ((*npoints) != last_npoints) {
      printf("interpolating at R = %e, eta = %e, delta = %e...\n",
	     p->R,p->eta,p->delta);
      copy_arrays(*r,*y,r_cp,y_cp,last_npoints); // copy arrays r and y into r_cp and y_cp
      free_matrices(c,s,y,r,rf_,integrand1,integrand2, // resize all arrays to new npoints size
		    last_npoints,ns);
      allocate_matrices(c,s,y,r,rf_,integrand1,integrand2,
			(*npoints),ns);
      propagate_r(*r,h,(*npoints));
      interpolate_array(*r,*y,r_cp,y_cp,(*npoints)); // interpolate old y values (now stored in y_cp)
      // so that the new y array has an initial guess for solvde.
      printf("done interpolating.\n");
      
    }
    
    // if last loop was successful with last_xpoint sized arrays, then just
    // use last y values as initial guess
    else propagate_r(*r,h,(*npoints));
    
    successful_solvde = solvde(itmax,conv,slowc,scalv,
			       ns,(*npoints),*y,*r,*c,*s,p,h); // relax to compute psi,
      //                                                            psi' curves, 
      


    if (!successful_solvde) {
      printf("solvde convergence failed, trying one more time with a "
	     "linear guess and a final twist angle value of pi/4");
      linearGuess(*r,*y,slopeguess,h,(*npoints));
      successful_solvde = solvde(itmax,conv,slowc,scalv,
				 ns,(*npoints),*y,*r,*c,*s,p,h); // relax to compute psi,
      //                                                            psi' curves, 
    }
    if (!successful_solvde) {
      write_failure("SOLVDE",*r,*y,*rf_,*npoints,*p);
    return;
    }
    
    // calculate energy, derivatives (see energy.c for code)
    successful_qromb = energy_stuff(E,dEdx,p,*r,*y,*rf_,
				    *integrand1,*integrand2,(*npoints));
    
    
    (*npoints) = ((*npoints)-1)*2+1;
    
  }
  while (!successful_qromb && (*npoints) <= max_size);

  (*npoints) = ((*npoints)-1)/2+1;



  if (!successful_qromb) { // writes the psi(r), r*f, etc, and exits
    write_failure("QROMB",*r,*y,*rf_,*npoints,*p);
  }

  return;
}
  

  
double backtracker(double beta,double rate,double E,double *dEdx,
		   struct params *p,double *r,double **y,
		   double *r_cp,double **y_cp,double conv,
		   int itmax,int npoints,struct arr_ns *ns,
		   int last_npoints,int max_size,double min_rate)
// Implements the standard backtracking routine for gradient descent. If         //
// x = (R,eta,delta)', E(x) is the energy, and grad = (ddR,ddeta,dddelta)', then //
// searching for rate such that E(x-rate*grad(E(x)))>E(x)-rate/2*grad(E(x))^2 is //
// satisfied. This algorithm gives me both the appropriate rate, and then the    //
// next values of R, eta, and delta. //
{

  double *rc;
  double **yc,**sc;
  double ***cc;
  double *rf_c, *integrand1c,*integrand2c;
  double Eold;
  double *tmpdEdx;
  double gradEsq;

  tmpdEdx = vector(1,3); 

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_arrays(r,y,rc,yc,last_npoints);

  // compute the current value of the gradient of E squared
  gradEsq = vector_norm(dEdx,3);

  // save the current value of E in Eold
  Eold = E;

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.

  update_p(p,rate,dEdx);
  
  // calculate E, given the new values of R,eta, and delta

  // we don't care about the new derivatives (they don't come into the algorithm),
  // so can just use temp variables to calculate them



  array_constant(0,tmpdEdx,3);
  // setting derivatives to zero flags single_calc function to not waste time
  // calculating them.

  single_calc(&E,tmpdEdx,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
	      &integrand2c,y_cp,r_cp,conv,itmax,&npoints,
	      last_npoints,ns,max_size);

  last_npoints = npoints;



  reset_p(p,rate,dEdx);

  while (E>=Eold-0.5*rate*gradEsq && rate > min_rate) {

    rate *=beta;

    update_p(p,rate,dEdx);
    single_calc(&E,tmpdEdx,p,&cc,&sc,&yc,&rc,&rf_c,&integrand1c,
		&integrand2c,y_cp,r_cp,conv,itmax,&npoints,
		last_npoints,ns,max_size);

    last_npoints = npoints;

    reset_p(p,rate,dEdx);

  }

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);
  free_vector(tmpdEdx,1,3);
  return rate;
}
/*
double beta_k()
{
  double beta;
  beta = (dEdR-lastdEdR)*dEdR+(dEdeta-lastdEdeta)*dEdeta+(dEddelta-lastdEddelta)*dEddelta;
  beta = dEdR*dEdR+dEdeta*dEdeta+dEddelta*dEddelta;


}


void set_p_k()
{


}
*/

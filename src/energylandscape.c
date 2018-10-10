#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"



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
			double **dEdvarlast,char scan_what[],double *R, 
			double *dEdR,double *dEdRlast,double *eta,
			double *dEdeta, double *dEdetalast,double *delta,
			double *dEddelta,double *dEddeltalast);

void allocate_matrices(double ****c,double ***s,double ***y,double **r,
		       double **rf_, double **integrand1,
		       double **integrand2,int npoints,struct arr_ns *ns);
void free_matrices(double ****c,double ***s,double ***y,double **r,
		   double **rf_, double **integrand1,
		   double **integrand2,int npoints,struct arr_ns *ns);
void single_calc(double *E,double *dEdR,double *dEdeta,double *dEddelta,
		 struct params *p,double ****c,double ***s,double ***y,
		 double **r,double **rf_,double **integrand1,
		 double **integrand2,double **y_cp,double *r_cp,
		 double conv,int itmax,int *npoints,int last_npoints,
		 struct arr_ns *ns,int max_size);
double backtracker(double beta,double rate,double E,double dEdR,double dEdeta,
		   double dEddelta,struct params *p,double *r,double **y,
		   double *r_cp, double **y_cp,double conv, int itmax,
		   int npoints,struct arr_ns *ns,int last_npoints,int max_size,
		   double min_rate);



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
			double **dEdvarlast,char scan_what[],double *R, 
			double *dEdR,double *dEdRlast,double *eta,
			double *dEdeta, double *dEdetalast,double *delta,
			double *dEddelta,double *dEddeltalast)
// convenient way to make a variables (ending in var) that have the  //
// same addresses as whichever parameter that is specified by        //
// scan_what.
{

  if (strcmp(scan_what,"R")==0) {
    printf("R!\n");
    *var = R;
    *var0 = *R;
    *dEdvar = dEdR;
    *dEdvarlast = dEdRlast;
  }
  else if (strcmp(scan_what,"eta")==0) {
    printf("eta!\n");
    *var = eta;
    *var0 = *eta;
    *dEdvar = dEdeta;
    *dEdvarlast = dEdetalast;
  }
  else if (strcmp(scan_what,"delta")==0) {
    printf("delta!\n");
    *var = delta;
    *var0 = *delta;
    *dEdvar = dEddelta;
    *dEdvarlast = dEddeltalast;
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
  double initialSlope;
  double E;
  double dEdR, dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar, *dEdvarlast;
  double Emin = 1e100;
  int max_size = (mpt-1)*4+1;
  double dx;
  int count_x;
  struct arr_ns ns;

  assign_ns(&ns);

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);


  setup_var_pointers(&var,&var0,&dEdvar,&dEdvarlast,scan_what,&p.R, 
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,
		     &p.delta,&dEddelta,&dEddeltalast);
  dx = (p.upperbound_x-var0)/num_x;
  printf("dx = %lf\n",dx);
  
  h = p.R/(npoints-1);
  // initial guess for the functional form of psi(r) and psi'(r)
  initialSlope = M_PI/(4.0*p.R);
  linearGuess(r,y,initialSlope,h,npoints); //linear initial guess 

  count_x = 0;

  while (count_x <= num_x) {

    // for each value of x (i.e *var) in E vs x


    single_calc(&E,&dEdR,&dEdeta,&dEddelta,&p,&c,&s,&y,&r,&rf_,&integrand1,
		&integrand2,y_cp,r_cp,conv,itmax,&npoints,last_npoints,
		&ns,max_size);

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
    dEdRlast = dEdR;
    dEdetalast = dEdeta;
    dEddeltalast = dEddelta;
    
  }

  free_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		npoints,&ns);
  free_matrix(y_cp,1,ns.nyj,1,max_size);
  free_vector(r_cp,1,max_size);

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
  double dEdR, dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar_x, *dEdvar_xlast;
  double *dEdvar_y, *dEdvar_ylast;
  double Emin = 1e100;
  int max_size = (mpt-1)*10+1;
  double dx,dy;
  int count_y,count_x;
  struct arr_ns ns;

  assign_ns(&ns);

  // initialize the pointers to the x variable (in E vs x vs y) so that
  // they reference the correct derivatives of E. x is either R,
  // eta, or delta.

  // malloc the relevant arrays
  allocate_matrices(&c,&s,&y,&r,&rf_,&integrand1,&integrand2,
		    npoints,&ns);
  // as well as two extra arrays to copy r and y contents into
  y_cp = matrix(1,ns.nyj,1,max_size);
  r_cp = vector(1,max_size);
  

  setup_var_pointers(&var_x,&var_x0,&dEdvar_x,&dEdvar_xlast,scan_what_x,&p.R, 
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,&p.delta,
		     &dEddelta,&dEddeltalast);
  dx = (p.upperbound_x-var_x0)/(num_x-1);

  setup_var_pointers(&var_y,&var_y0,&dEdvar_y,&dEdvar_ylast,scan_what_y,&p.R,
		     &dEdR,&dEdRlast,&p.eta,&dEdeta,&dEdetalast,&p.delta,
		     &dEddelta,&dEddeltalast);
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

      single_calc(&E,&dEdR,&dEdeta,&dEddelta,&p,&c,&s,&y,&r,&rf_,&integrand1,
		  &integrand2,y_cp,r_cp,conv,itmax,&npoints,last_npoints,
		  &ns,max_size);

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

      dEdRlast = dEdR;
      dEdetalast = dEdeta;
      dEddeltalast = dEddelta;
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
  
  return;
}

void graddesc(struct params p,FILE *energy,FILE *psi,
	      FILE *denergydR,FILE *denergydeta,
	      FILE *denergyddelta,FILE *surfacetwist,
	      FILE *energydensity,double conv,int itmax,
	      int mpt,double rate)
{
  int npoints = mpt;
  int last_npoints = mpt;
  double h;
  double **y,**y_cp,*r,*r_cp, **s, ***c;
  double *rf_,*integrand1,*integrand2;
  double initialSlope;
  double E;
  double beta = 0.7;
  double Elast = -1e100; // a really large, negative number to start with
  double dEdR, dEdeta,dEddelta;
  double lastdEdR,lastdEdeta,lastdEddelta;
  double min_rate = 1e3*conv;
  double last_rate = rate;
  double lastR = p.R;
  double lasteta = p.eta;
  double lastdelta = p.delta;
  int max_size = (mpt-1)*8+1;
  int count = 0;
  bool small_rate = false;

  struct arr_ns ns;
  assign_ns(&ns);
  
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

  dEdR = dEdeta = dEddelta = lastdEdR = lastdEdeta = lastdEddelta = 1e100;


  while (fabs(dEdR) > conv || fabs(dEdeta) > conv
	 || fabs(dEddelta) > conv) {

    single_calc(&E,&dEdR,&dEdeta,&dEddelta,&p,&c,&s,&y,&r,&rf_,&integrand1,
		&integrand2,y_cp,r_cp,conv,itmax,&npoints,last_npoints,
		&ns,max_size);

    last_npoints = npoints;

    fprintf(energy,"%d\t%.8e\n",count,E);
    fprintf(denergydR,"%.8e\t%.8e\n",p.R,dEdR);
    fprintf(denergydeta,"%.8e\t%.8e\n",p.eta,dEdeta);
    fprintf(denergyddelta,"%.8e\t%.8e\n",p.delta,dEddelta);
    fprintf(surfacetwist,"%.8e\t%.8e\n",p.R,y[1][mpt]);


    last_rate = backtracker(beta,rate,E,dEdR,dEdeta,dEddelta,&p,r,y,r_cp,y_cp,
			    conv,itmax,npoints,&ns,last_npoints,max_size,min_rate);

    count += 1;
    printf("count = %d\n",count);

    if (last_rate <= min_rate) {
      printf("last rate <= min_rate, where min_rate = %e and last "
	     " rate = %e.\n",min_rate,last_rate);
      if (fabs(lastR-p.R)<conv && fabs(lasteta-p.eta)<conv
	  && fabs(lastdelta-p.delta)<conv) {
	printf("changes in R, eta, and delta are smaller than %e.\n",conv);
	
	if (dEdR*lastdEdR <=0 && dEdeta*lastdEdeta <= 0
	    && dEddelta*lastdEddelta <= 0) {
	  printf("jumped over the minimimum (derivatives all changed signs).\n");
	  break;
	}
      }
    }

    lastR = p.R;
    lasteta = p.eta;
    lastdelta = p.delta;
    lastdEdR = dEdR;
    lastdEdeta = dEdeta;
    lastdEddelta = dEddelta;
    Elast = E;

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
  
  return;
}

void single_calc(double *E,double *dEdR,double *dEdeta,double *dEddelta,
		 struct params *p,double ****c,double ***s,double ***y,
		 double **r,double **rf_,double **integrand1,
		 double **integrand2,double **y_cp,double *r_cp,
		 double conv,int itmax,int *npoints,int last_npoints,
		 struct arr_ns *ns,int max_size)
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
    successful_qromb = energy_stuff(E,dEdR,dEdeta,dEddelta,p,*r,*y,
				   *rf_,*integrand1,*integrand2,(*npoints));
    
    
    (*npoints) = ((*npoints)-1)*2+1;
    
  }
  while (!successful_qromb && (*npoints) <= max_size);

  (*npoints) = ((*npoints)-1)/2+1;



  if (!successful_qromb) { // writes the psi(r), r*f, etc, and exits
    write_failure("QROMB",*r,*y,*rf_,*npoints,*p);
  }

  return;
}
  
  
double backtracker(double beta,double rate,double E,double dEdR,double dEdeta,
		   double dEddelta,struct params *p,double *r,double **y,
		   double *r_cp, double **y_cp,double conv, int itmax,
		   int npoints,struct arr_ns *ns,int last_npoints,int max_size,
		   double min_rate)
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
  double tmpdEdR,tmpdEdeta,tmpdEddelta;
  double gradEsq;
  double rate0=rate;

  // create arrays to manipulate without affecting external arrays
  allocate_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  // copy the external r,y arrays into the internal rc,yc arrays
  copy_arrays(r,y,rc,yc,last_npoints);

  // compute the current value of the gradient of E squared
  gradEsq = dEdR*dEdR+dEdeta*dEdeta+dEddelta*dEddelta;

  // save the current value of E in Eold
  Eold = E;

  // calculate the next values of R, eta, and delta, if the rate was just the normal
  // rate input to this function.
  p->R = p->R-rate*dEdR;
  p->eta = p->eta-rate*dEdeta;
  p->delta = p->delta-rate*dEddelta;    
  
  // calculate E, given the new values of R,eta, and delta

  // we don't care about the new derivatives (they don't come into the algorithm),
  // so can just use temp variables to calculate them


  tmpdEdR = tmpdEdeta = tmpdEddelta = 0;

  // setting derivatives to zero flags single_calc function to not waste time
  // calculating them.


  single_calc(&E,&tmpdEdR,&tmpdEdeta,&tmpdEddelta,p,&cc,&sc,&yc,&rc,
	      &rf_c,&integrand1c,&integrand2c,y_cp,r_cp,conv,itmax,
	      &npoints,last_npoints,ns,max_size);

  last_npoints = npoints;


  p->R = p->R+rate*dEdR;
  p->eta = p->eta+rate*dEdeta;
  p->delta = p->delta+rate*dEddelta;    


  while (E>=Eold-0.5*rate*gradEsq && rate > min_rate) {

    rate *=beta;

    p->R = p->R-rate*dEdR;
    p->eta = p->eta-rate*dEdeta;
    p->delta = p->delta-rate*dEddelta;    

    single_calc(&E,&tmpdEdR,&tmpdEdeta,&tmpdEddelta,p,&cc,&sc,&yc,&rc,
		&rf_c,&integrand1c,&integrand2c,y_cp,r_cp,conv,itmax,
		&npoints,last_npoints,ns,max_size);

    last_npoints = npoints;

    p->R = p->R+rate*dEdR;
    p->eta = p->eta+rate*dEdeta;
    p->delta = p->delta+rate*dEddelta;


  }

  p->R = p->R-rate*dEdR;
  p->eta = p->eta-rate*dEdeta;
  p->delta = p->delta-rate*dEddelta;    

  if (rate0 != rate) printf("rate = %e, rate0 = %e\n",rate,rate0);

  free_matrices(&cc,&sc,&yc,&rc,&rf_c,&integrand1c,&integrand2c,npoints,ns);

  return rate;
}



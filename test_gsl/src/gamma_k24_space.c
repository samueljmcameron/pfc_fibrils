#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "../../shared_src/nrutil.h"
#include "headerfile.h"


int main(int argc, char **argv)
{
  void initialize_params(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_x(double *x,struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void linearGuess(double *r, double **y, double initialSlope,double h,int mpt);

  bool drive(double *E,struct params *p,double *x,FILE *energy);

  void set_scanning_values(double *gammalow,double *gammahigh,int *num_gamma,
			   double *k24low,double *k24high,int *num_k24,char **args);

  double gamma_j(int j, double gammalow, double gammahigh,int num_gamma);

  double k24_i(int i, double k24low, double k24high,int num_k24);
  
  void set_scalings(double Rtemp, double etatemp, double deltatemp,struct params *p);

  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);
  p.mpt = (MAX_M-1)/8+1;


  double *x;
  x = vector(1,X_SIZE);
  initialize_x(x,&p);

  //make a linear initial guess for psi(r)
  double h = x[1]/(p.mpt-1);
  double slopeguess = M_PI/(4.0*x[1]);
  linearGuess(p.r,p.y,slopeguess,h,p.mpt);


  FILE *energy;
  initialize_file(&energy,argv[1],"energy",p);

  FILE *surfacetwist;
  initialize_file(&surfacetwist,argv[1],"surfacetwist",p);

  FILE *radius;
  initialize_file(&radius,argv[1],"radius",p);

  FILE *eta;
  initialize_file(&eta,argv[1],"eta",p);

  FILE *delta;
  initialize_file(&delta,argv[1],"delta",p);

  double gammalow,gammahigh,k24low,k24high;
  int num_gamma,num_k24;
  set_scanning_values(&gammalow,&gammahigh,&num_gamma,
		      &k24low,&k24high,&num_k24,argv);

  // temp_i variables are stored values from i-1 iteration. Since k24
  // is the i iterate, this means everytime gamma = gammalow, the guess
  // variable used is from k24 = (k24high-k24low)/(num_k24-1)*(i-1).
  double Rtemp_i = x[1];
  double etatemp_i = x[2];
  double deltatemp_i = x[3];


  // temp_j variables are stored values from j-1 iteration. Since gamma
  // is the j iterate, this means everytime the scan is moving at constant
  // k24, the guess variable is used from the previous gamma's solution.
  double Rtemp_j = x[1];
  double etatemp_j = x[2];
  double deltatemp_j = x[3];





  double E;
  int i,j;

  for (i = 0; i < num_k24; i++) {

    p.k24 = k24_i(i,k24low,k24high,num_k24);


    if (i != 0) {

      // if not first calculation, use xvalues at the 0,i-1 grid point, as well
      // as a linear guess for the form of y[1][1..p.mpt] = psi(r) and
      // y[2][1..p.mpt] = psi'(r).

      set_scalings(Rtemp_i,etatemp_i,deltatemp_i,&p);
      p.mpt = (MAX_M-1)/8+1;
      h = Rtemp_i/(p.mpt-1);
      linearGuess(p.r,p.y,slopeguess,h,p.mpt);
    }

    for (j = 0; j < num_gamma; j++) {

      p.gamma_s = gamma_j(j,gammalow,gammahigh,num_gamma);


      // if not first calculation at the current i value, use xvalues at 
      // the j-1,i grid point
      if (j != 0) set_scalings(Rtemp_j,etatemp_j,deltatemp_j,&p);
      
      initialize_x(x,&p);
  
      if (drive(&E,&p,x,(NULL))) {
	
	if (j == 0) {

	  // if at gamma=gammalow, then save the equilibrium x values, and
	  // the initial slope of psi(r), to use as guesses for 0,i+1.

	  Rtemp_i = x[1];
	  etatemp_i = x[2];
	  deltatemp_i = x[3];
	  slopeguess = p.y[2][1];
	}
	
	
	Rtemp_j = x[1];
	etatemp_j = x[2];
	deltatemp_j = x[3];

	fprintf(energy,"%13.6e\t",E);
	fprintf(surfacetwist,"%13.6e\t",p.y[1][p.mpt]);
	fprintf(radius,"%13.6e\t",x[1]);
	fprintf(eta,"%13.6e\t",x[2]);
	fprintf(delta,"%13.6e\t",x[3]);

      } else {

	fprintf(energy,"%13.6e\t",sqrt(-1));
	fprintf(surfacetwist,"%13.6e\t",sqrt(-1));
	fprintf(radius,"%13.6e\t",sqrt(-1));
	fprintf(eta,"%13.6e\t",sqrt(-1));
	fprintf(delta,"%13.6e\t",sqrt(-1));

      }

    }

    fprintf(energy,"\n");
    fprintf(surfacetwist,"\n");
    fprintf(radius,"\n");
    fprintf(eta,"\n");
    fprintf(delta,"\n");

  }

  free_vector(x,1,X_SIZE);
  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(energy);
  fclose(surfacetwist);
  fclose(radius);
  fclose(eta);
  fclose(delta);

  return 0;

}

void set_scalings(double Rtemp, double etatemp, double deltatemp,struct params *p)
{
  p->Rguess = Rtemp;
  p->Rupper = Rtemp*1.5;
  p->Rlower = Rtemp*0.75;
  p->etaguess = etatemp;
  p->etaupper = etatemp+0.1;
  p->etalower = etatemp-0.02;
  p->deltaguess = deltatemp;
  p->deltaupper = 0.817;
  p->deltalower = 0.95*deltatemp;

  return;
}

double gamma_j(int j, double gammalow, double gammahigh,int num_gamma)
{
  double dgamma = (gammahigh-gammalow)/(num_gamma-1);
  return j*dgamma+gammalow;
}


double k24_i(int i, double k24low, double k24high,int num_k24)
{
  double dk24 = (k24high-k24low)/(num_k24-1);
  return i*dk24+k24low;
}

void set_scanning_values(double *gammalow,double *gammahigh,int *num_gamma,
			 double *k24low,double *k24high,int *num_k24,char **args)
{
  sscanf(args[17],"%lf",gammalow);
  sscanf(args[18],"%lf",gammahigh);
  sscanf(args[19],"%d",num_gamma);
  sscanf(args[20],"%lf",k24low);
  sscanf(args[21],"%lf",k24high);
  sscanf(args[22],"%d",num_k24);
  return;
}


void initialize_x(double *x,struct params *p)
{
  x[1] = p->Rguess;
  x[2] = p->etaguess;
  x[3] = p->deltaguess;
  return;
}

void initialize_param_vectors(struct params *p)
{

  p->r = vector(1,MAX_M);
  p->y = matrix(1,NE,1,MAX_M);
  p->r_cp = vector(1,MAX_M);
  p->y_cp = matrix(1,NE,1,MAX_M);
  p->rf_fib = vector(1,MAX_M);
  p->s = matrix(1,NSI,1,NSJ);
  p->c = f3tensor(1,NCI,1,NCJ,1,MAX_M+1);

  return;
}

void initialize_params(struct params *p,char **args)
{


  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->d0);
  sscanf(args[6],"%lf",&p->omega);
  sscanf(args[7],"%lf",&p->gamma_s);
  sscanf(args[8],"%lf",&p->Rguess);
  sscanf(args[9],"%lf",&p->etaguess);
  sscanf(args[10],"%lf",&p->deltaguess);
  sscanf(args[11],"%lf",&p->Rupper);
  sscanf(args[12],"%lf",&p->Rlower);
  sscanf(args[13],"%lf",&p->etaupper);
  sscanf(args[14],"%lf",&p->etalower);
  sscanf(args[15],"%lf",&p->deltaupper);
  sscanf(args[16],"%lf",&p->deltalower);

  printf("parameter values for minimization:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
  printf("d0 = %e\n",p->d0);
  printf("omega = %e\n",p->omega);
  printf("gamma_s = %e\n",p->gamma_s);


  printf("initial guesses for R, eta, and delta:\n");

  printf("guess for R = %e\n",p->Rguess);
  printf("guess for eta = %e\n",p->etaguess);
  printf("guess for delta = %e\n",p->deltaguess);


  printf("estimated ranges of R, eta, and delta:\n");

  printf("upper bound estimate of R = %e\n",p->Rupper);
  printf("lower bound estimate of R = %e\n",p->Rlower);
  printf("upper bound estimate of eta = %e\n",p->etaupper);
  printf("lower bound estimate of eta = %e\n",p->etalower);
  printf("upper bound estimate of delta = %e\n",p->deltaupper);
  printf("lower bound estimate of delta = %e\n",p->deltalower);

  return;

}

void initialize_file(FILE **output,char *path,char *fname,struct params p)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_"
	   "%1.4e.txt",p.K33,p.Lambda,p.d0,p.omega);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}

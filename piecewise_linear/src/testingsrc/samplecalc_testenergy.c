#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "funcs/nrutil.h"
#include "headerfile.h"


int main(int argc, char **argv)
{

  void setparams(struct params *p,char **args);

  double Efunc(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void save_psivsr(FILE *psivsr,double *r,double **y,double *rf_fib,int n);

  void propagate_r(double *r, double h,int mpt);

  void buildpsi(double *r,double **y,struct params *p,int n);

  void compute_rf2233b1(double *rf_fib,struct params *p,double *r, double **y, int n);

  struct params p;

  int n = 2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2+1;


  double *r,*rf_fib;
  double **y;

  r = vector(1,n);
  rf_fib = vector(1,n);

  y = matrix(1,2,1,n);
  
  setparams(&p,argv);

  double h = p.R/(n-1);

  printf("E discrete (C code) from discrete integral = %.8lf\n",Efunc(&p));

  FILE *psi_r;

  initialize_file(&psi_r,argv[1],"psivsr_discrete",p);

  propagate_r(r,h,n);

  buildpsi(r,y,&p,n);

  compute_rf2233b1(rf_fib,&p,r,y,n);

  save_psivsr(psi_r,r,y,rf_fib,n);

  fclose(psi_r);

  free_vector(r,1,n);
  free_vector(rf_fib,1,n);
  free_matrix(y,1,2,1,n);
  
  return 0;

}

void save_psivsr(FILE *psivsr,double *r,double **y,double *rf_fib,int n)
{
  int i;

  for (i = 1; i <= n; i++) {
    fprintf(psivsr,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",r[i],y[1][i],y[2][i],
	    rf_fib[i]);
  }
  return;
}


void setparams(struct params *p,char **args)
{
  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->R);
  sscanf(args[8],"%lf",&p->eta);
  sscanf(args[9],"%lf",&p->delta);
  sscanf(args[10],"%lf",&p->R_c);
  sscanf(args[11],"%lf",&p->R_s);
  sscanf(args[12],"%lf",&p->psip_c);
  sscanf(args[13],"%lf",&p->psip_s);
  sscanf(args[14],"%lf",&p->psip_R);
  
  return;
}

void initialize_file(FILE **output,char *path,char *fname,struct params p)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e.txt",p.K33,p.k24,p.Lambda,p.omega,p.gamma_s);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}

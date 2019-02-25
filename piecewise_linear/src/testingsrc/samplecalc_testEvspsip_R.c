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

  double dEdpsip_R(struct params *p);

  double ufunc(double x_1,double x_2,double zeta);
  double dudx_1(double x_1,double zeta);
  double dudx_2(double x_2,double zeta);
  double dudzeta(double x_1,double x_2,double zeta);
  
  double vfunc(double x_1,double x_2,double xi,double zeta);
  double dvdx_1(double x_1,double xi,double zeta);
  double dvdx_2(double x_2,double xi,double zeta);
  double dvdxi(double x_1,double x_2,double xi, double zeta);
  double dvdzeta(double x_1,double x_2,double xi,double zeta);
  
  double f_1func(double x_1,double x_2,double xi,double zeta);
  double df_1dx_1(double x_1,double xi,double zeta);
  double df_1dx_2(double x_2,double xi,double zeta);
  double df_1dxi(double x_1,double x_2,double xi,double zeta);
  double df_1dzeta(double x_1,double x_2,double xi,double zeta);
  
  double f_2func(double x_1,double x_2,double xi,double zeta);
  double df_2dx_1(double x_1,double xi,double zeta);
  double df_2dx_2(double x_2,double xi,double zeta);
  double df_2dxi(double x_1,double x_2,double xi,double zeta);
  double df_2dzeta(double x_1,double x_2,double xi,double zeta);
  
  double g_1func(double x_1,double x_2,double xi,double zeta);
  double dg_1dx_1(double x_1,double xi,double zeta);
  double dg_1dx_2(double x_2,double xi,double zeta);
  double dg_1dxi(double x_1,double x_2,double xi,double zeta);
  double dg_1dzeta(double x_1,double x_2,double xi,double zeta);
  
  double g_2func(double x_1,double x_2,double xi,double zeta);
  double dg_2dx_1(double x_1,double xi,double zeta);
  double dg_2dx_2(double x_2,double xi,double zeta);
  double dg_2dxi(double x_1,double x_2,double xi,double zeta);
  double dg_2dzeta(double x_1,double x_2,double xi,double zeta);
  
  
  void initialize_file(FILE **output,char *path,char *fname,struct params p);


  struct params p;

  setparams(&p,argv);
  
  FILE *Evspsip_R;
  
  initialize_file(&Evspsip_R,argv[1],"Evspsip_R",p);

  double psip_R0 = -0.0006;
  double psip_Rf = 0.0006;
  int num = 60;

  double dpsip_R = (psip_Rf-psip_R0)/num;
  
  int i;

  double psi2;

  for (i = 0; i < num; i++) {

    p.psip_R = dpsip_R*i+psip_R0;

    fprintf(Evspsip_R,"%13.6e\t%13.6e\t%13.6e\t",p.psip_R,Efunc(&p),dEdpsip_R(&p));

    psi2 = (p.psip_s-p.psip_R)*p.R_s+(p.psip_c-p.psip_s)*p.R_c;

    fprintf(Evspsip_R,"%13.12lf\n",dg_1dzeta(p.R_s,p.R,psi2,p.psip_R));

  }

  fclose(Evspsip_R);
  
  return 0;

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

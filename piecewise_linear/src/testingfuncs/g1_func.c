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

  
  double g_1func(double x_1,double x_2,double xi,double zeta);
  double dg_1dx_1(double x_1,double xi,double zeta);
  double dg_1dx_2(double x_2,double xi,double zeta);
  double dg_1dxi(double x_1,double x_2,double xi,double zeta);
  double dg_1dzeta(double x_1,double x_2,double xi,double zeta);
    
  void initialize_file(FILE **output,char *path,char *fname,struct params p);


  struct params p;

  setparams(&p,argv);
  
  FILE *g1_func;
  
  initialize_file(&g1_func,argv[1],"g1_func",p);

  double psip_R0 = -1;
  double psip_Rf = 1;
  int num = 200;

  double dpsip_R = (psip_Rf-psip_R0)/num;

  double x_1,x_2,xi,zeta;
  
  int i;

  double psi2;

  for (i = 0; i < num; i++) {

    p.psip_R = dpsip_R*i+psip_R0;

    x_1 = p.R_s;
    x_2 = p.R;
    zeta = p.psip_R;
    xi = (p.psip_s-p.psip_R)*p.R_s+(p.psip_c-p.psip_s)*p.R_c;

    fprintf(g1_func,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",zeta,
	    g_1func(x_1,x_2,xi,zeta),dg_1dxi(x_1,x_2,xi,zeta),
	    dg_1dzeta(x_1,x_2,xi,zeta));


  }

  fclose(g1_func);
  
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

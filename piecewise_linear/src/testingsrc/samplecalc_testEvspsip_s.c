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

  double dEdpsip_s(struct params *p);
  
  void initialize_file(FILE **output,char *path,char *fname,struct params p);


  struct params p;

  setparams(&p,argv);
  
  FILE *Evspsip_s;
  
  initialize_file(&Evspsip_s,argv[1],"Evspsip_s",p);

  double psip_s0 = -0.05;
  double psip_sf = .95;
  int num = 1000;

  double dpsip_s = (psip_sf-psip_s0)/num;
  
  int i;

  for (i = 0; i < 1000; i++) {

    p.psip_s = dpsip_s*i+psip_s0;

    fprintf(Evspsip_s,"%13.6e\t%13.6e\t%13.6e\n",p.psip_s,Efunc(&p),dEdpsip_s(&p));

  }

  fclose(Evspsip_s);
  
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

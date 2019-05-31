#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "energy_src/nrutil.h"
#include "headerfile.h"

int main(int argc, char **argv)
{
  void initialize_params(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  double E_calc(struct params *p);

  void initialize_first_xpoint_x0(double *x0,struct params *p);

  void initialize_second_xpoint_x1(double *x1,struct params *p);

  void x_of_t(struct params *p,const double t,const double *x0,const double *x1);

  double set_t_start(char **args);

  double set_t_end(char **args);

  int set_numpoints(char **args);

  void initialize_psi_and_r(struct params *p);

  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  double *x0; // this will hold the coordinates at the first minimum (smaller R)
  p.x_size = 3;
  x0 = vector(1,p.x_size);
  initialize_first_xpoint_x0(x0,&p);

  double *x1; // this will hold the coordinates at the second minimum (larger R)
  x1 = vector(1,p.x_size);
  initialize_second_xpoint_x1(x1,&p);


  FILE *Evst;
  initialize_file(&Evst,argv[1],"Evst",p);
  
  double t_start = set_t_start(argv);
  double t_end = set_t_end(argv);
  int numpoints = set_numpoints(argv);
  double t_step = (t_end-t_start)/(numpoints-1);

  double t;
  double E;


  int iter;

  x_of_t(&p,t_start,x0,x1);

  initialize_psi_and_r(&p);


  for (iter = 0; iter < numpoints; iter++) {

    t = t_start + iter*t_step;

    x_of_t(&p,t,x0,x1);

    E = E_calc(&p);
    printf("E=%lf\n",E);

    fprintf(Evst,"%13.6e\t%13.9e\n",t,E);

  }

  free_vector(x0,1,p.x_size);
  free_vector(x1,1,p.x_size);
  
  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(Evst);

  return 0;

}

int set_numpoints(char **args)
{
  int n;
  sscanf(args[15],"%d",&n);
  return n;
}

double set_t_start(char **args)
{
  double t0;
  sscanf(args[13],"%lf",&t0);
  return t0;
}

double set_t_end(char **args)
{
  double t1;
  sscanf(args[14],"%lf",&t1);
  return t1;
}


void x_of_t(struct params *p,const double t,const double *x0,const double *x1)
{
  int i;

  p->R = x0[1]+t*(x1[1]-x0[1]);
  p->eta = x0[2]+t*(x1[2]-x0[2]);
  p->delta = x0[3]+t*(x1[3]-x0[3]);
  
  
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

void initialize_psi_and_r(struct params *p)
{

  void linearGuess(double *r, double **y, double initialSlope,double h,int mpt);

  double h = p->R/(p->mpt-1);
  double slopeguess = M_PI/(4.0*p->R);
  linearGuess(p->r,p->y,slopeguess,h,p->mpt); //linear initial guess for psi(r)
  
  return;
}
  
void initialize_params(struct params *p,char **args)
{


  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->R0);
  sscanf(args[8],"%lf",&p->R1);
  sscanf(args[9],"%lf",&p->eta0);
  sscanf(args[10],"%lf",&p->eta1);
  sscanf(args[11],"%lf",&p->delta0);
  sscanf(args[12],"%lf",&p->delta1);
  p->mpt = (MAX_M-1)/8+1;

  printf("parameter values for calculation:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
  printf("omega = %e\n",p->omega);
  printf("gamma_s = %e\n",p->gamma_s);

  printf("the two x points along the line x(t)=x0+t*v are:\n");

  printf("x0 = (%e,%e,%e)'\n",p->R0,p->eta0,p->delta0);
  printf("x1 = (%e,%e,%e)'\n",p->R1,p->eta1,p->delta1);
  return;

}

void initialize_first_xpoint_x0(double *x0,struct params *p)
{
  x0[1] = p->R0;
  x0[2] = p->eta0;
  x0[3] = p->delta0;
  return;
}

void initialize_second_xpoint_x1(double *x1,struct params *p)
{
  x1[1] = p->R1;
  x1[2] = p->eta1;
  x1[3] = p->delta1;
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

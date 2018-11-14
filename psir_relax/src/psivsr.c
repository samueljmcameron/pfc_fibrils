#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../../src/nrutil.h"
#include "../../src/headerfile.h"

#define NE 2                   // # of 1st order DEs
#define M 2*2*2*2*2*2*2*2*2*2*2*2*2*2+1  // # of mesh points (2^M+1 for romberg integration)
#define NB 1                   // # of BCs at first boundary (k = 1)
#define NSI NE                 // max # i of S_i,j
#define NSJ (2*NE+1)           // max # j of S_i,j
#define NYJ NE                 // # of dependent variables
#define NYK M                  // number of points in each dependent variable
#define NCI NE                 // # of rows in storage matrix within c[m][:][:]
#define NCJ (NE-NB+1)          // # of columns in storage matrix within c[m][:][:]
#define NCK (M+1)              // # number of points in tensor c, c[:][m][n]

#define X_SIZE     3
#define ITMAX      10000
#define CONV_ODE   1e-10
#define CONV_MIN   1e-8
#define MAX_SIZE_M (8)*(M-1)+(1)
#define EPS        1e-14

int main(int argc, char **argv)
{
  struct params p;
  char scan_what[20];
  char path[200];
  char suffix[200];
  char f1[200];
  FILE *psivsr;
  double psip01,psip02;
  double *x;

  x = vector(1,X_SIZE);

  // read in the the variables and paths
  snprintf(path,sizeof(path),"%s",argv[1]);
  sscanf(argv[2],"%lf",&p.K33);
  sscanf(argv[3],"%lf",&p.k24);
  sscanf(argv[4],"%lf",&p.Lambda);
  sscanf(argv[5],"%lf",&p.d0);
  sscanf(argv[6],"%lf",&p.omega);
  sscanf(argv[7],"%lf",&x[1]);
  sscanf(argv[8],"%lf",&x[2]);
  sscanf(argv[9],"%lf",&x[3]);
  sscanf(argv[10],"%lf",&p.gamma_s);


  printf("K33 = %lf\n",p.K33);
  printf("k24 = %lf\n",p.k24);
  printf("Lambda = %lf\n",p.Lambda);
  printf("d0 = %lf\n",p.d0);
  printf("omega = %lf\n",p.omega);
  printf("R = %lf\n",x[1]);
  printf("eta = %lf\n",x[2]);
  printf("delta = %lf\n",x[3]);
  printf("gamma_s = %lf\n",p.gamma_s);
  

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,x[1],x[2],x[3],p.gamma_s);

  snprintf(f1,sizeof(f1),"%s_psivsr_%s",
	   path,suffix);

  psivsr = fopen(f1,"w");

  double *r;
  double **y;
  double psip0;
  double h;

  
  int itmax = 1000;
  double conv = 1e-14;
  double *scalv;
  struct arr_ns ns;
  double ***c;
  double **s;
  double *r_cp,*rf_fib;
  double **y_cp;
  double E;
  bool solv = false;
  bool Ecalc;

  assign_ns(&ns);

  scalv = vector(1,2);
  scalv[1] = 0.1;
  scalv[2] = 4.0;

  allocate_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,M);

  h = x[1]/(M-1);

  propagate_r(r,h,M);

  psip01 = 0.0;
  psip02 = M_PI/(0.01*x[1]);
  psip0 = brent(psip01,psip02,EPS,1000,r,y,&p,x,h,M);

  printf("before solvde, y[2][1] = %.12e, y[1][M] = %e\n",y[2][1],y[1][M]);
  printf("f = %.12e\n",F_eq(y[1][M],y[2][M],x[1],&p));
  Ecalc = successful_E_count(&E,&p,x,r,y,rf_fib,M);
  printf("E calc success: %d, E = %e\n",Ecalc,E);
  solv = solvde(itmax,conv,scalv,&ns,M,r,y,c,s,&p,x,h);
  printf("after solvde, y[2][1] = %.12e, y[1][M] = %e\n",y[2][1],y[1][M]);
  printf("f = %.12e\n",F_eq(y[1][M],y[2][M],x[1],&p));
  Ecalc = successful_E_count(&E,&p,x,r,y,rf_fib,M);
  printf("E calc success: %d, E = %e\n",Ecalc,E);
  if (!solv) {
    printf("calculation failed!\n");
    exit(1);
  }
  
  psip0 = y[2][1];
  printf("psip0 =  %e\n",psip0);

  save_psi(psivsr,r,y,M);

  free_matrices(ns,&c,&s,&r,&y,&r_cp,&y_cp,&rf_fib,M);

  
  free_vector(x,1,X_SIZE);

  fclose(psivsr); // close file!

  return 0;
}

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../../src/nrutil.h"
#include "../../src/headerfile.h"

#define NE 2                   // # of 1st order DEs
#define M 2*2*2*2*2*2*2*2*2*2*2*2+1  // # of mesh points (2^M+1 for romberg integration)
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

int main(int argc, char **argv)
{
  struct params p;
  char scan_what[20];
  char path[200];
  char suffix[200];
  char f1[200];
  FILE *bc;
  double psip01,psip02;
  int numpoints;
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
  sscanf(argv[11],"%lf",&psip01);
  sscanf(argv[12],"%lf",&psip02);
  sscanf(argv[13],"%d",&numpoints);


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

  snprintf(f1,sizeof(f1),"%s_bcvspsip0_%s",
	   path,suffix);

  bc = fopen(f1,"w");

  shootscan_driver(p,x,bc,psip01,psip02,numpoints,M);
  

  free_vector(x,1,X_SIZE);

  fclose(bc); // close file!

  return 0;
}

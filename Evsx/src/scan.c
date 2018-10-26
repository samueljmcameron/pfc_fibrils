#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../../src/nrutil.h"
#include "../../src/headerfile.h"

#define NE 2                   // # of 1st order DEs
#define M 2*2*2*2*2*2*2*2*2*2*2+1  // # of mesh points (2^M+1 for romberg integration)
#define NB 1                   // # of BCs at first boundary (k = 1)
#define NSI NE                 // max # i of S_i,j
#define NSJ (2*NE+1)           // max # j of S_i,j
#define NYJ NE                 // # of dependent variables
#define NYK M                  // number of points in each dependent variable
#define NCI NE                 // # of rows in storage matrix within c[m][:][:]
#define NCJ (NE-NB+1)          // # of columns in storage matrix within c[m][:][:]
#define NCK (M+1)              // # number of points in tensor c, c[:][m][n]

#define X_SIZE 3               // size of vector x = (R,eta,delta)



int main(int argc, char **argv)
{
  int itmax = 10000;
  double conv = 1.0e-10;
  struct params p;
  char scan_what[20];
  char path[200];
  char suffix[200];
  char f1[200],f2[200];
  FILE *energy, *psi;
  int num_scan = 200;
  double *x;
  int x_index;
  x = vector(1,X_SIZE); // x = (R,eta,delta)'

  
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
  sscanf(argv[11],"%lf",&p.upperbound_x);

  printf("K33 = %lf\n",p.K33);
  printf("k24 = %lf\n",p.k24);
  printf("Lambda = %lf\n",p.Lambda);
  printf("d0 = %lf\n",p.d0);
  printf("omega = %lf\n",p.omega);
  printf("R = %lf\n",x[1]);
  printf("eta = %lf\n",x[2]);
  printf("delta = %lf\n",x[3]);
  printf("gamma_s = %lf\n",p.gamma_s);
  printf("upperbound_x = %lf\n",p.upperbound_x);
  printf("upperbound_y = %lf\n",p.upperbound_y);
  
  snprintf(scan_what,sizeof(scan_what),"%s",argv[12]);

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,x[1],x[2],x[3],p.gamma_s,
	   p.upperbound_x);

  snprintf(f1,sizeof(f1),"%s_%s_energy_%s",
	   path,scan_what,suffix);
  snprintf(f2,sizeof(f2),"%s_%s_psivsr_%s",
	   path,scan_what,suffix);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");

  x_index = index(scan_what);
  
  scanE(p,x,energy,psi,conv,itmax,M,num_x,x_index);
  
  free_vector(x,1,X_SIZE);
  fclose(energy); // close file!
  fclose(psi);

  return 0;
}

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

int main(int argc, char **argv)
{


  double **y,*r, **s, ***c;

  int itmax= 1000,logStep = 1000;
  double conv = 1.0e-10;
  double K33;
  double k24;
  double Lambda;
  double d0;
  double omega;
  double R;
  double L;
  double eta;
  double delta;
  double gamma_s;
  double gamma_t;
  double upperbound;
  double initialSlope;
  char scan_what[20];
  char path[200];
  char f1[200],f2[200];
  FILE *energy, *psi;

  // read in the the variables and paths
  snprintf(path,sizeof(path),"%s",argv[1]);
  sscanf(argv[2],"%lf",&K33);
  sscanf(argv[3],"%lf",&k24);
  sscanf(argv[4],"%lf",&Lambda);
  sscanf(argv[5],"%lf",&d0);
  sscanf(argv[6],"%lf",&omega);
  sscanf(argv[7],"%lf",&R);
  sscanf(argv[8],"%lf",&L);
  sscanf(argv[9],"%lf",&eta);
  sscanf(argv[10],"%lf",&delta);
  sscanf(argv[11],"%lf",&gamma_s);
  sscanf(argv[12],"%lf",&gamma_t);
  sscanf(argv[13],"%lf",&upperbound);
  snprintf(scan_what,sizeof(scan_what),"%s",argv[14]);

  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);
  r = vector(1,NYK);

  initialSlope = M_PI/(4.0*R);

  snprintf(f1,sizeof(f1),"%s_%s_energy_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   path,scan_what,K33,k24,Lambda,d0,omega,R,L,eta,delta,
	   gamma_s,gamma_t,upperbound);
  snprintf(f2,sizeof(f2),"%s_%s_psivsr_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   path,scan_what,K33,k24,Lambda,d0,omega,R,L,eta,delta,
	   gamma_s,gamma_t,upperbound);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");


  scanE(r,y,c,s,K33,k24,Lambda,d0,omega,R,L,eta,delta,gamma_s,
	gamma_t,initialSlope,energy,psi,conv,itmax,M,upperbound,
	scan_what);
  

  fclose(energy); // close file!
  fclose(psi);
  free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
  free_matrix(s,1,NSI,1,NSJ);
  free_matrix(y,1,NYJ,1,NYK);
  free_vector(r,1,NYK);
  return 0;
}

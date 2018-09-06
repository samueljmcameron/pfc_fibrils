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
  int itmax = 10000;
  double conv = 1.0e-10;
  struct params p;
  double K33;
  double k24;
  double Lambda;
  double d0;
  double omega;
  double R;
  double eta;
  double delta;
  double gamma_s;
  double upperbound;
  char scan_what[20];
  char path[200];
  char suffix[200];
  char f1[200],f2[200];
  FILE *energy, *psi;

  // read in the the variables and paths
  snprintf(path,sizeof(path),"%s",argv[1]);
  sscanf(argv[2],"%lf",&p.K33);
  sscanf(argv[3],"%lf",&p.k24);
  sscanf(argv[4],"%lf",&p.Lambda);
  sscanf(argv[5],"%lf",&p.d0);
  sscanf(argv[6],"%lf",&p.omega);
  sscanf(argv[7],"%lf",&p.R);
  sscanf(argv[8],"%lf",&p.eta);
  sscanf(argv[9],"%lf",&p.delta);
  sscanf(argv[10],"%lf",&p.gamma_s);
  sscanf(argv[11],"%lf",&p.upperbound);
  K33 = p.K33;
  k24 = p.k24;
  Lambda = p.Lambda;
  d0 = p.d0;
  omega = p.omega;
  R = p.R;
  eta = p.eta;
  delta = p.delta;
  gamma_s = p.gamma_s;
  upperbound = p.upperbound;
  /*
  sscanf(argv[2],"%lf",&K33);
  sscanf(argv[3],"%lf",&k24);
  sscanf(argv[4],"%lf",&Lambda);
  sscanf(argv[5],"%lf",&d0);
  sscanf(argv[6],"%lf",&omega);
  sscanf(argv[7],"%lf",&R);
  sscanf(argv[8],"%lf",&eta);
  sscanf(argv[9],"%lf",&delta);
  sscanf(argv[10],"%lf",&gamma_s);
  sscanf(argv[11],"%lf",&upperbound);
  */
  snprintf(scan_what,sizeof(scan_what),"%s",argv[12]);

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,p.delta,p.gamma_s,
	   p.upperbound);

  snprintf(f1,sizeof(f1),"%s_%s_energy_%s",
	   path,scan_what,suffix);
  snprintf(f2,sizeof(f2),"%s_%s_psivsr_%s",
	   path,scan_what,suffix);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");

  scanE(p,energy,psi,conv,itmax,M,scan_what);

  //  scanE(K33,k24,Lambda,d0,omega,R,eta,delta,gamma_s,
  //energy,psi,conv,itmax,M,upperbound,scan_what);
  

  fclose(energy); // close file!
  fclose(psi);

  return 0;
}

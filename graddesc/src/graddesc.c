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
  char scan_what[20];
  char path[200];
  char suffix[200];
  char f1[200],f2[200],f3[200],f4[200],f5[200],f6[200];
  FILE *energy,*psi,*denergydR,*denergydeta;
  FILE *denergyddelta,*surfacetwist;
  double rateR, rateeta, ratedelta;

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
  sscanf(argv[11],"%lf",&rateR);
  sscanf(argv[12],"%lf",&rateeta);
  sscanf(argv[13],"%lf",&ratedelta);

  printf("K33 = %lf\n",p.K33);
  printf("k24 = %lf\n",p.k24);
  printf("Lambda = %lf\n",p.Lambda);
  printf("d0 = %lf\n",p.d0);
  printf("omega = %lf\n",p.omega);
  printf("R = %lf\n",p.R);
  printf("eta = %lf\n",p.eta);
  printf("delta = %lf\n",p.delta);
  printf("gamma_s = %lf\n",p.gamma_s);
  printf("rateR = %lf\n",rateR);
  printf("rateeta = %lf\n",rateeta);
  printf("ratedelta = %lf\n",ratedelta);
  
  snprintf(scan_what,sizeof(scan_what),"%s",argv[12]);

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,p.delta,p.gamma_s);

  snprintf(f1,sizeof(f1),"%s_energy_%s",
	   path,suffix);
  snprintf(f2,sizeof(f2),"%s_psivsr_%s",
	   path,suffix);
  snprintf(f3,sizeof(f3),"%s_dEdR_%s",
	   path,suffix);
  snprintf(f4,sizeof(f4),"%s_dEdeta_%s",
	   path,suffix);
  snprintf(f5,sizeof(f5),"%s_dEddelta_%s",
	   path,suffix);
  snprintf(f6,sizeof(f6),"%s_surfacetwist_%s",
	   path,suffix);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");
  denergydR = fopen(f3,"w");
  denergydeta = fopen(f4,"w");
  denergyddelta = fopen(f5,"w");
  surfacetwist = fopen(f6,"w");

  graddesc(p,energy,psi,denergydR,denergydeta,denergyddelta,
	   surfacetwist,conv,itmax,M,rateR,rateeta,ratedelta);
  

  fclose(energy); // close file!
  fclose(psi);
  fclose(denergydR);
  fclose(denergydeta);
  fclose(denergyddelta);
  fclose(surfacetwist);

  return 0;
}
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

  int itmax= 10000;
  double conv = 1.0e-10;
  struct params p;
  char scan_what_x[20],scan_what_y[20];
  char path[200];
  char suffix[200],f1[200],f2[200],f3[200];
  char f4[200],f5[200];
  FILE *energy, *psi, *deriv_energy_x;
  FILE *deriv_energy_y, *surfacetwist;
  int num_x = 200, num_y = 201;

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
  sscanf(argv[11],"%lf",&p.upperbound_x);
  sscanf(argv[12],"%lf",&p.upperbound_y);
  snprintf(scan_what_x,sizeof(scan_what_x),"%s",argv[13]);
  snprintf(scan_what_y,sizeof(scan_what_y),"%s",argv[14]);

  printf("K33 = %lf\n",p.K33);
  printf("k24 = %lf\n",p.k24);
  printf("Lambda = %lf\n",p.Lambda);
  printf("d0 = %lf\n",p.d0);
  printf("omega = %lf\n",p.omega);
  printf("R = %lf\n",p.R);
  printf("eta = %lf\n",p.eta);
  printf("delta = %lf\n",p.delta);
  printf("gamma_s = %lf\n",p.gamma_s);
  printf("upperbound_x = %lf\n",p.upperbound_x);
  printf("upperbound_y = %lf\n",p.upperbound_y);

  printf("scan_what_x = %s\n",scan_what_x);
  printf("scan_what_y = %s\n",scan_what_y);


  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   p.K33,p.k24,p.Lambda,p.d0,p.omega,p.R,p.eta,p.delta,p.gamma_s,
	   p.upperbound_x,p.upperbound_y);
  
  snprintf(f1,sizeof(f1),"%s_%s_%s_energy_%s",
	   path,scan_what_x,scan_what_y,suffix);

  snprintf(f2,sizeof(f2),"%s_%s_%s_psivsr_%s",
	   path,scan_what_x,scan_what_y,suffix);
  
  snprintf(f3,sizeof(f3),"%s_%s_%s_deriv_energy_%s_%s",
	   path,scan_what_x,scan_what_y,scan_what_x,suffix);

  snprintf(f4,sizeof(f4),"%s_%s_%s_deriv_energy_%s_%s",
	   path,scan_what_x,scan_what_y,scan_what_y,suffix);

  snprintf(f5,sizeof(f5),"%s_%s_%s_surfacetwist_%s",
	   path,scan_what_x,scan_what_y,suffix);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");
  deriv_energy_x = fopen(f3,"w");
  deriv_energy_y = fopen(f4,"w");
  surfacetwist = fopen(f5,"w");
  


  scan2dE(p,energy,psi,deriv_energy_x,deriv_energy_y,
	  surfacetwist,conv,itmax,M,num_x,num_y,
	  scan_what_x,scan_what_y);
  
  

  fclose(energy); // close file!
  fclose(psi);
  fclose(deriv_energy_x);
  fclose(deriv_energy_y);
  fclose(surfacetwist);
  return 0;
}

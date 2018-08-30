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
  double eta;
  double delta;
  double gamma_s;
  double upperbound_x,upperbound_y;
  char scan_what_x[20],scan_what_y[20];
  char path[200];
  char suffix[200],f1[200],f2[200],f3[200];
  char f4[200],f5[200];
  FILE *energy, *psi, *deriv_energy_x;
  FILE *deriv_energy_y, *surfacetwist;

  // read in the the variables and paths
  snprintf(path,sizeof(path),"%s",argv[1]);
  sscanf(argv[2],"%lf",&K33);
  sscanf(argv[3],"%lf",&k24);
  sscanf(argv[4],"%lf",&Lambda);
  sscanf(argv[5],"%lf",&d0);
  sscanf(argv[6],"%lf",&omega);
  sscanf(argv[7],"%lf",&R);
  sscanf(argv[8],"%lf",&eta);
  sscanf(argv[9],"%lf",&delta);
  sscanf(argv[10],"%lf",&gamma_s);
  sscanf(argv[11],"%lf",&upperbound_x);
  sscanf(argv[12],"%lf",&upperbound_y);
  snprintf(scan_what_x,sizeof(scan_what_x),"%s",argv[13]);
  snprintf(scan_what_y,sizeof(scan_what_y),"%s",argv[14]);

  printf("K33 = %lf\n",K33);
  printf("k24 = %lf\n",k24);
  printf("Lambda = %lf\n",Lambda);
  printf("d0 = %lf\n",d0);
  printf("omega = %lf\n",omega);
  printf("R = %lf\n",R);
  printf("eta = %lf\n",eta);
  printf("delta = %lf\n",delta);
  printf("gamma_s = %lf\n",gamma_s);
  printf("upperbound_x = %lf\n",upperbound_x);
  printf("upperbound_y = %lf\n",upperbound_y);

  printf("scan_what_x = %s\n",scan_what_x);
  printf("scan_what_y = %s\n",scan_what_y);



  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);
  r = vector(1,NYK);



  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e_%1.4e_%1.4e_%1.4e_%1.4e.txt",
	   K33,k24,Lambda,d0,omega,R,eta,delta,gamma_s,
	   upperbound_x,upperbound_y);
  
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
  


  scan2dE(r,y,c,s,K33,k24,Lambda,d0,omega,R,eta,delta,gamma_s,
	  energy,psi,deriv_energy_x,deriv_energy_y,
	  surfacetwist,conv,itmax,M,upperbound_x,upperbound_y,
	  scan_what_x,scan_what_y);
  
  

  fclose(energy); // close file!
  fclose(psi);
  fclose(deriv_energy_x);
  fclose(deriv_energy_y);
  fclose(surfacetwist);
  free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
  free_matrix(s,1,NSI,1,NSJ);
  free_matrix(y,1,NYJ,1,NYK);
  free_vector(r,1,NYK);
  return 0;
}

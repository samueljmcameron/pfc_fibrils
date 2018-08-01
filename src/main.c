#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"
#include "headerfile.h"

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
  void solvde(int itmax, double conv, double slowc, double scalv[],
	      int ne, int nb, int m, double **y, double *r,
	      double ***c, double **s, double K33, double k24,
	      double Lambda, double eta, double d0, double L,
	      double h);

  double **y,*r, **s, ***c;
  double K33 = 30.0;
  double k24 = 0.8;
  double gamma_s = 0.05;
  double gamma_t = .01;

  int itmax= 1000,logStep = 1000;
  double conv = 1.0e-10;
  double Lambda = 10.0;
  double d0;
  double R;
  double L;
  double eta;
  char scan_what[20];

  double initialSlope;

  char path[200];
  char f1[200],f2[200];
  FILE *energy, *psi;

  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);
  r = vector(1,NYK);

  sscanf(argv[1],"%lf",&d0);
  sscanf(argv[2],"%lf",&R);
  sscanf(argv[3],"%lf",&L);
  sscanf(argv[4],"%lf",&eta);

  snprintf(scan_what,sizeof(scan_what),"%s",argv[5]);

  K33 = 30.;
  initialSlope = M_PI/(4.0*R);

  snprintf(path,sizeof(path),"../data/");
  snprintf(f1,sizeof(f1),"%sd0_%1.4e_R_%1.4e_eta_%1.4e_Evs%s.txt",
	   path,d0,R,eta,scan_what);
  snprintf(f2,sizeof(f2),"%sd0_%1.4e_R_%1.4e_eta_%1.4e_psiVSr_%s.txt",
	   path,d0,R,eta,scan_what);

  energy = fopen(f1,"w");
  psi = fopen(f2,"w");

  scanE(r,y,c,s,K33,k24,Lambda,eta,d0,L,R,initialSlope,
	 gamma_s,gamma_t,energy,psi,conv,itmax,M,scan_what);


  fclose(energy); // close file!
  fclose(psi);
  free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
  free_matrix(s,1,NSI,1,NSJ);
  free_matrix(y,1,NYJ,1,NYK);
  free_vector(r,1,NYK);
  return 0;
}

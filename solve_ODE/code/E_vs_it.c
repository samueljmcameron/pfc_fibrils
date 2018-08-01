#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"

#define NE 2                   // # of 1st order DEs
#define M 2*2*2*2*2*2*2*2*2+1  // # of mesh points (2^M+1 for romberg integration)
#define NB 1                   // # of BCs at first boundary (k = 1)
#define NSI NE                 // max # i of S_i,j
#define NSJ (2*NE+1)           // max # j of S_i,j
#define NYJ NE                 // # of dependent variables
#define NYK M                  // number of points in each dependent variable
#define NCI NE                 // # of rows in storage matrix within c[m][:][:]
#define NCJ (NE-NB+1)          // # of columns in storage matrix within c[m][:][:]
#define NCK (M+1)              // # number of points in tensor c, c[:][m][n]


void linearGuess(double *r, double **y, double initialSlope, double h);

int main(int argc, char **argv)
{
  void solvde(int itmax, double conv, double slowc, double scalv[],
	      int ne, int nb, int m, double **y, double *r,
	      double ***c, double **s, double K33, double k24,
	      double Lambda, double eta, double d0, double L,
	      double h, int mpt);

  double **y,*r, **s, ***c;
  double h,K33,k24;

  int itmax= 1000;
  double conv = 1.0e-10;
  double eta = 1.0;
  double Lambda;
  double d0 = 1.0;
  double L = 1.0;
  double R;

  double initialSlope;

  double slowc = 1.0;
  double scalv[NE+1];
  
  scalv[1] = .1;    // guess for the order of magnitude of the psi values
  scalv[2] = 4.0;   // guess for the order of magnitude of the psi' values



  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);
  r = vector(1,NYK);

  k24 = 0.8;
  sscanf(argv[1],"%lf",&Lambda);
  sscanf(argv[2],"%lf",&R);

  K33 = 30.;

  initialSlope = 1e-4*M_PI/(4.0*R);
  printf("initialSlope = %e\n",initialSlope);

  h = R/(M-1);

  linearGuess(r,y,initialSlope,h); // if the first curve being computed, 
  //                                  use linear initial guess
    
  solvde(itmax,conv,slowc,scalv,NE,NB,M,y,r,c,s,K33,k24,
	 Lambda,eta,d0,L,h,M); // relax to compute psi while plotting
  //                              E vs iteration


  free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
  free_matrix(s,1,NSI,1,NSJ);
  free_matrix(y,1,NYJ,1,NYK);
  free_vector(r,1,NYK);
  return 0;
}





void linearGuess(double *r, double **y, double initialSlope, double h)
{
  int k;
  
  for (k=1;k <=M; k++) { // initial guess!
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*r[k]; // y1 is psi
    y[2][k] = initialSlope; // y2 is psi'!!!!!!!!!!!!!!!!!!
  }
  return;
}

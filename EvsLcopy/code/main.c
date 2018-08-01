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

void linearGuess(double *r, double **y, double initialSlope, double h);

void propagate_r(double *r, double h);

void savePsi(FILE *Psi,double *r, double **y);

void saveEnergy(FILE *f, double R, double Energy, double derivative,
		double observable);

void grad_desc(double *r,double **y,double ***c,double **s,
	       double K33,double k24,double Lambda,
	       double eta,double d0,double L,double R,
	       double initialSlope,double gamma_s,
	       double gamma_t,FILE *f,FILE *Psi,double conv,
	       int itmax,char search_what[]);

int main(int argc, char **argv)
{
  void solvde(int itmax, double conv, double slowc, double scalv[],
	      int ne, int nb, int m, double **y, double *r,
	      double ***c, double **s, double K33, double k24,
	      double Lambda, double eta, double d0, double L,
	      double h, int mpt);

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
  char search_what[20];

  double initialSlope;

  char path[200];
  char f1[200],f2[200];
  FILE *f, *Psi;

  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);
  r = vector(1,NYK);

  sscanf(argv[1],"%lf",&d0);
  sscanf(argv[2],"%lf",&R);
  sscanf(argv[3],"%lf",&L);
  sscanf(argv[4],"%lf",&eta);

  snprintf(search_what,sizeof(search_what),"%s",argv[5]);

  K33 = 30.;
  initialSlope = M_PI/(4.0*R);

  snprintf(path,sizeof(path),"../Results/");
  snprintf(f1,sizeof(f1),"%sd0_%1.4e_R_%1.4e_eta_%1.4e_Evs%s.txt",
	   path,d0,R,eta,search_what);
  snprintf(f2,sizeof(f2),"%sd0_%1.4e_R_%1.4e_eta_%1.4e_PsiVSr_%s.txt",
	   path,d0,R,eta,search_what);

  f = fopen(f1,"w");
  Psi = fopen(f2,"w");

  grad_desc(r,y,c,s,K33,k24,Lambda,eta,d0,L,R,initialSlope,
	    gamma_s,gamma_t,f,Psi,conv,itmax,search_what);


  fclose(f); // close file!
  fclose(Psi);
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

void propagate_r(double *r, double h)
{
  int k;
  for (k=1;k <=M; k++) r[k] = (k-1)*h; // only change r since psi, psi' are stored from last loop
  return;
}

void savePsi(FILE *Psi,double *r, double **y)
{
  int i;

  for (i = 1; i <= M; i++) {
    fprintf(Psi,"%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i]);
  }
  printf("psi(R) = %1.2e\n",y[1][M]);
  return;
}

void saveEnergy(FILE *f, double R, double Energy, double derivative,
		double observable)
{
  fprintf(f,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,Energy,
	  derivative,observable);
  return;
}


void grad_desc(double *r,double **y,double ***c,double **s,
	       double K33,double k24,double Lambda,
	       double eta,double d0,double L,double R,
	       double initialSlope,double gamma_s,
	       double gamma_t,FILE *f,FILE *Psi,double conv,
	       int itmax,char search_what[])
{
  int isitone=1;
  double *var,var0;
  double upperbound;
  double h;
  double slowc = 1.0;
  double scalv[NE+1];
  double rf_[M+1],integrand1[M+1],integrand2[M+1];
  double E;
  double dEdR, dEdL, dEdeta;
  double dEdRlast = dEdR;
  double dEdLlast = dEdL;
  double dEdetalast = dEdeta;
  double *dEdvar, *dEdvarlast;

  if (strcmp(search_what,"R")==0) {
    printf("R!\n");
    var = &R;
    var0 = R;
    upperbound = 3.0;
    dEdvar = &dEdR;
    dEdvarlast = &dEdRlast;
  }
  else if (strcmp(search_what,"L")==0) {
    printf("L!\n");
    var = &L;
    var0 = L;
    upperbound = 40.0;
    dEdvar = &dEdL;
    dEdvarlast = &dEdLlast;
  }
  else if (strcmp(search_what,"eta")==0) {
    printf("eta!\n");
    var = &eta;
    var0 = eta;
    dEdvar = &dEdeta;
    dEdvarlast = &dEdetalast;
  }
  else {
    printf("Need either R, L, or eta as argv[5] input."
	   "Exiting to system.\n");
    exit(1);
  }

  scalv[1] = .1;    // guess for the order of magnitude of the psi values
  scalv[2] = 4.0;   // guess for the order of magnitude of the psi' values
  
  while (*var <= upperbound) {

    h = R/(M-1);

    if (isitone == 1) {
      linearGuess(r,y,initialSlope,h); // use linear initial guess 
      isitone += 1;
    }
    
    else propagate_r(r,h); // if not initial value (not first loop), 
    //                          use previous result as initial guess

    solvde(itmax,conv,slowc,scalv,NE,NB,M,y,r,c,s,K33,k24,
	   Lambda,eta,d0,L,h,M); // relax to compute psi,
    // psi' curves

    energy_stuff(&E,&dEdR,&dEdeta,&dEdL,R,K33,k24,Lambda,eta,d0,L,gamma_s,
		 gamma_t,r,y,rf_,integrand1,integrand2,M);
    
    saveEnergy(f,*var,E,*dEdvar,y[1][M]);   // save L,E, and surface twist

    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0 && *dEdvarlast < 0) {
      savePsi(Psi,r,y);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E+0.5); 
    }

    *var += 0.001;
    dEdRlast = dEdR;
    dEdLlast = dEdL;
    dEdetalast = dEdeta;
    
  }

  return;
}

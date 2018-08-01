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

void compute_ffibril23(double K33, double K24, double *r,
		       double **y, double *rf_);

double E1_plus_E2(double R,double *r,double *rf_);

double E_R(double R, double K24, double gamma,double *r,
	   double **y, double *rf_);

double dEdR(double R, double K33, double K24,
	    double *r,double **y, double *rf_,double gamma);

void linearGuess(double *r, double **y, double initialSlope, double h);

void propagate_r(double *r, double h);

void savePsi(FILE *Psi,double *r, double **y);

void saveEnergy(FILE *f, double R, double Energy, double derivative,
		double twistprimeR);

void grad_desc_EvsR(double *r, double *rf_, double **y, double ***c,double **s,
		    double initialSlope,double gamma,FILE *f, FILE *Psi,
		    double conv, int itmax);

void search_EvsR(double *r, double *rf_, double **y, double ***c,double **s,
		  double initialSlope,double gamma,FILE *f, FILE *Psi,
		  double expRlast, double conv, int itmax);

int mpt=M;
double h, r[M+1],K33,K24;  // global variables communicating with difeq


int main(int argc, char **argv)
{
  void solvde(int itmax, double conv, double slowc, double scalv[],
	      int indexv[], int ne, int nb, int m, double **y, double ***c, double **s);

  double **y, **s, ***c;
  double rf_[M+1];

  int itmax= 1000,logStep = 1000;

  double conv = 1.0e-8;

  double initialSlope,expRlast;

  double gamma;

  char path[200];
  char f1[200],f2[200];
  FILE *f, *Psi;

  y = matrix(1,NYJ,1,NYK);
  s = matrix(1,NSI,1,NSJ);
  c = f3tensor(1,NCI,1,NCJ,1,NCK);

  expRlast = 0.5;

  sscanf(argv[1],"%lf",&K24);
  sscanf(argv[2],"%lf",&gamma);

  K33 = 30.;

  initialSlope = 100;

  snprintf(path,sizeof(path),"../Results/");
  snprintf(f1,sizeof(f1),"%sK24is_%1.4e_GAMMAis%1.4eEvsR.txt",
	     path,K24,gamma);
  snprintf(f2,sizeof(f2),"%sK24is_%1.4e_GAMMAis%1.4ePsiVSr.txt",
	   path,K24,gamma);

  f = fopen(f1,"w");
  Psi = fopen(f2,"w");

  // now going to loop over range of R values, R = 10^expR with expR varying as necessary

  grad_desc_EvsR(r,rf_,y,c,s,initialSlope,gamma,f,Psi,conv,itmax);
  //  search_EvsR(r,rf_,y,c,s,initialSlope,gamma,f,Psi,expRlast,conv,itmax);

  fclose(f); // close file!
  fclose(Psi);
  free_f3tensor(c,1,NCI,1,NCJ,1,NCK);
  free_matrix(s,1,NSI,1,NSJ);
  free_matrix(y,1,NYJ,1,NYK);
  return 0;
}

void compute_ffibril23(double K33, double K24, double *r,
		       double **y, double *rf_)
{
  int i;
  double siny, sin2y;
  
  rf_[1] = 0;
  
  for (i = 2; i <= M; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    rf_[i] = 0.5*(r[i]+y[2][i]*y[2][i]*r[i]+0.5*0.5*sin2y	\
		      *sin2y/r[i]-2*y[2][i]*r[i]-sin2y+y[2][i]*sin2y) \
      +0.5*K33*(siny*siny*siny*siny)/r[i];//-sin2y*K24*y[2][i];
  }
  return;
}

double E1_plus_E2(double R,double *r,double *rf_)
{
  return 2/(R*R)*qromb(r,rf_,M);
}

double E_R(double R, double K24, double gamma,double *r,
	   double **y, double *rf_)
{

  double E;
  E = E1_plus_E2(R,r,rf_); // first calculate bulk energy per unit length
  E = E+2.0/R*(-K24*(sin(y[1][M])*sin(y[1][M]))/R+gamma);  // add surface term tension term

  return E;
}

double dEdR(double R, double K33, double K24,
	    double *r,double **y, double *rf_,double gamma)
{
  double ans;
  double sin2y = sin(2*y[1][M]);
  double siny = sin(y[1][M]);

  ans = -2.0/R*E1_plus_E2(R,r,rf_);

  ans += 2.0/R*(0.5*(1-y[2][M]-sin2y/(2*R))*
		(1-y[2][M]-sin2y/(2*R))+
		0.5*K33*(siny*siny*siny*siny)/(R*R));
  ans += 2.0/(R*R)*(2.0*K24*(siny*siny/R-siny*cos(y[1][M])*y[2][M])-gamma);

  return ans;
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
		double twistprimeR)
{
  fprintf(f,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,Energy,
	  derivative,twistprimeR); // store (E,R) value in file
  return;
}

void grad_desc_EvsR(double *r, double *rf_, double **y, double ***c,double **s,
		    double initialSlope,double gamma,FILE *f, FILE *Psi,
		    double conv, int itmax)
{
  int isitone=1, firsttime=1;
  int indexv[NE+1];
  double R;
  double rate = 0;
  double Rdot=1e10;
  double slowc = 1.0;
  double scalv[NE+1];
  double lastEnergy = 1e15;
  double Energy;
  double derivative = 1e15;
  
  indexv[1] = 1;
  indexv[2] = 2;
  scalv[1] = .1;    // guess for the order of magnitude of the psi values
  scalv[2] = 4.0;   // guess for the order of magnitude of the psi' values

  R = pow(10,-3);

  while (fabs(derivative) > 1e-8) {
  
    h = R/(M-1);

    if (isitone == 1) {
      linearGuess(r,y,initialSlope,h); // use linear initial guess 
      isitone += 1;
    }
    
    else propagate_r(r,h); // if expR is not initial value (not first loop), 
    //                          use previous result as initial guess

    solvde(itmax,conv,slowc,scalv,indexv,NE,NB,M,y,c,s); // relax to compute psi,
                                                         // psi' curves    
    compute_ffibril23(K33,K24,r,y,rf_);  // compute r*f_fibril, store in rf_ array

    Energy = E_R(R,K24,gamma,r,y,rf_);

    derivative = dEdR(R,K33,K24,r,y,rf_,gamma);

    if (rate == 0) {
      rate = fabs(0.1/derivative);
    }
    else rate *= 1.1;
    Rdot = -rate*derivative;
    
    saveEnergy(f,R,Energy,derivative,y[2][M]);   // save R,E, and surface twist

    R = Rdot+R;


  }


  printf("%e\t%e\n",rate,R);
  printf("%e\n",derivative);
  printf("y[2][M]/y[2][1] = %lf\n",y[2][M]/y[2][1]);
  printf("R = %lf\n",R);
  
  savePsi(Psi,r,y);
  printf("SAVED!\n");
  printf("Energy at minimum = %1.2e\n",Energy);
  

  return;
}

void search_EvsR(double *r, double *rf_, double **y, double ***c,double **s,
		 double initialSlope,double gamma,FILE *f, FILE *Psi,
		 double expRlast, double conv, int itmax)
{
  int i;
  int isitone, firsttime=1;
  int indexv[NE+1];
  double R;
  double expR;
  double slowc = 1.0;
  double scalv[NE+1];
  double lastEnergy = 1e15;
  double Energy;
  double derivative;
  
  indexv[1] = 1;
  indexv[2] = 2;
  scalv[1] = .1;    // guess for the order of magnitude of the psi values
  scalv[2] = 4.0;   // guess for the order of magnitude of the psi' values
  

  for (isitone = 1, expR=-3; expR<=expRlast;expR += .0005,isitone++) { // log steps for plt
    R = pow(10,expR);
    h = R/(M-1);

    if (isitone == 1) linearGuess(r,y,initialSlope,h); // if the first curve being computed, 
    //                                                         use linear initial guess
    
    else propagate_r(r,h); // if expR is not initial value (not first loop), 
    //                          use previous result as initial guess

    //    printf("%lf\n",R);
    solvde(itmax,conv,slowc,scalv,indexv,NE,NB,M,y,c,s); // relax to compute psi,

    compute_ffibril23(K33,K24,r,y,rf_);  // compute r*f_fibril, store in rf_ array

    Energy = E_R(R,K24,gamma,r,y,rf_);

    derivative = dEdR(R,K33,K24,
		      r,y,rf_,gamma);

    saveEnergy(f,R,Energy,derivative,y[2][M]);   // save R,E, and surface twist

    if (Energy > lastEnergy && R > 1e-4 && firsttime == 1) {

      firsttime = 0;
      printf("y[2][M]/y[2][1] = %lf\n",y[2][M]/y[2][1]);
      printf("R = %lf\n",R);

      //if (fabs(y[2][M]/y[2][1]-10)<0.005) {
      savePsi(Psi,r,y);
      printf("SAVED!\n");
	// }
      printf("Energy at minimum = %1.2e\n",Energy);

    }
    lastEnergy = Energy;
  }
  //  savePsi(Psi,r,y);
  return;
}

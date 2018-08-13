#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"
#include "headerfile.h"

void linearGuess(double *r, double **y, double initialSlope,
		 double h,int mpt)
{
  int k;
  
  for (k=1;k <=mpt; k++) { // initial guess!
    r[k] = (k-1)*h;
    y[1][k] = initialSlope*r[k]; // y1 is psi
    y[2][k] = initialSlope; // y2 is psi'!!!!!!!!!!!!!!!!!!
  }
  return;
}

void propagate_r(double *r, double h,int mpt)
{
  int k;
  for (k=1;k <=mpt; k++) r[k] = (k-1)*h; // only change r since psi, psi' are stored from last loop
  return;
}

void save_psi(FILE *psi,double *r, double **y,int mpt)
{
  int i;

  for (i = 1; i <= mpt; i++) {
    fprintf(psi,"%10.8e\t%10.8e\t%10.8e\n",r[i],y[1][i],y[2][i]);
  }
  printf("psi(R) = %1.2e\n",y[1][mpt]);
  return;
}

void saveEnergy(FILE *energy, double R, double E, double derivative,
		double observable)
{
  fprintf(energy,"%10.8e\t%10.8e\t%10.8e\t%10.8e\n",R,E,
	  derivative,observable);
  return;
}


void scanE(double *r,double **y,double ***c,double **s,
	   double K33,double k24,double Lambda,double d0,
	   double omega,double R,double L,double eta,
	   double delta,double gamma_s,double gamma_t,
	   double initialSlope,FILE *energy,FILE *psi,
	   double conv,int itmax,int mpt, 
	   double upperbound, char scan_what[])
// The energy of the system E(R,L,eta). This function //
// generates data of E vs scan_what[] (either "L",    //
// "R", or "eta"), while holding the other two values //
// constant. The E vs scan_what[] data is saved in    //
// energy file. If a minimum (or multiple minima) are //
// found in the energy landscape, psi(r) is saved to  //
// the psi file. //
{
  int isitone=1;
  double *var,var0;
  double h;
  double slowc = 1.0;
  double scalv[2+1];
  double rf_[mpt+1],integrand1[mpt+1],integrand2[mpt+1];
  double E;
  double dEdR, dEdL, dEdeta,dEddelta;
  double dEdRlast = dEdR;
  double dEdLlast = dEdL;
  double dEdetalast = dEdeta;
  double dEddeltalast = dEddelta;
  double *dEdvar, *dEdvarlast;
  double Emin = 1e100;

  if (strcmp(scan_what,"R")==0) {
    printf("R!\n");
    var = &R;
    var0 = R;
    dEdvar = &dEdR;
    dEdvarlast = &dEdRlast;
  }
  else if (strcmp(scan_what,"L")==0) {
    printf("L!\n");
    var = &L;
    var0 = L;
    dEdvar = &dEdL;
    dEdvarlast = &dEdLlast;
  }
  else if (strcmp(scan_what,"eta")==0) {
    printf("eta!\n");
    var = &eta;
    var0 = eta;
    dEdvar = &dEdeta;
    dEdvarlast = &dEdetalast;
  }
  else if (strcmp(scan_what,"delta")==0) {
    printf("/delta!\n");
    var = &delta;
    var0 = delta;
    dEdvar = &dEddelta;
    dEdvarlast = &dEddeltalast;
  }
  else {
    printf("Need either R, L, eta, or delta as argv[14] input."
	   "Exiting to system.\n");
    exit(1);
  }

  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values
  
  while (*var <= upperbound) {

    h = R/(mpt-1);

    if (isitone == 1) {
      linearGuess(r,y,initialSlope,h,mpt); //linear initial guess 
      isitone += 1;
    }
    
    else propagate_r(r,h,mpt); // if not first loop, previous 
    //                            result is initial guess

    solvde(itmax,conv,slowc,scalv,2,1,mpt,y,r,c,s,K33,k24,
	   Lambda,d0,L,eta,delta,h); // relax to compute psi,
    //                                  psi' curves, note the 2,1
    //                                  corresponds to two eqns,
    //                                   and 1 BC at the r = 0.

    // calculate energy, derivatives (see energy.c for code)
    energy_stuff(&E,&dEdR,&dEdL,&dEdeta,&dEddelta,K33,k24,
		 Lambda,d0,omega,R,L,eta,delta,gamma_s,gamma_t,
		 r,y,rf_,integrand1,integrand2,mpt);


    // save L,E, and surface twist
    saveEnergy(energy,*var,E,*dEdvar,y[1][mpt]);

    if (*var != var0 && *dEdvar*(*dEdvarlast) <= 0
	&& *dEdvarlast < 0 && E <= Emin) {
      save_psi(psi,r,y,mpt);
      printf("SAVED!\n");
      printf("E_min-E_chol = %1.2e\n",E+0.5);
      Emin = E;
    }

    *var += 0.001;
    dEdRlast = dEdR;
    dEdLlast = dEdL;
    dEdetalast = dEdeta;
    dEddeltalast = dEddelta;
    
  }

  return;
}

#ifndef _HEADERFILE_H_
#define _HEADERFILE_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include <stdio.h>
#include <stdlib.h>


#define ZERO 1e-14
#define SMALL 1e-8
#define ABSERR 1e-10
#define RELERR 1e-10
#define SIZE_INTEGRATION 1000

#define FAILED_E 1e300


struct params{
  // these five parameters completely specify the parameter
  // space which corresponds to physical changes of fibril
  // environment.
  double K33;
  double k24;
  double Lambda;
  double omega;
  double gamma_s;

  // these next 24 parameters are just hyperparameters which 
  // only help with computation.

  // initial guesses for the gradient descent
  double R_guess;
  double eta_guess;
  double delta_guess;
  double R_c_guess;
  double R_s_guess;
  double psip_c_guess;
  double psip_s_guess;
  double psip_R_guess;

  // estimated upper and lower bounds for the gradient descent.
  double R_upper;
  double R_lower;
  double eta_upper;
  double eta_lower;
  double delta_upper;
  double delta_lower;
  double R_c_upper;
  double R_c_lower;
  double R_s_upper;
  double R_s_lower;
  double psip_c_upper;
  double psip_c_lower;
  double psip_s_upper;
  double psip_s_lower;
  double psip_R_upper;
  double psip_R_lower;
};

#ifndef _HEADERFILE_H_
#define _HEADERFILE_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include <stdio.h>
#include <stdlib.h>


#define ZERO 1e-14
#define SMALL 5e-3
#define ABSERR 1e-12
#define RELERR 1e-12
#define SIZE_INTEGRATION 1000

#define CONV_ODE   1e-10
#define CONV_MIN   1e-8

#define FAILED_E   1e300

#define DRIVER_FAILURE 0
#define DRIVER_SUCCESS 1
#define DRIVER_POORSCALING 2

#define DELTA_CLOSE_TO_ZERO 1e-5


struct params{
  // these five parameters completely specify the parameter
  // space which corresponds to physical changes of fibril
  // environment.
  double K33;
  double k24;
  double Lambda;
  double omega;
  double gamma_s;

  // these rest of parameters below are just hyperparameters which 
  // only help with computation.

  int x_size;

  double Escale;
  
  // initial guesses for the gradient descent
  double Rguess;
  double etaguess;
  double deltaguess;
  double R_cguess;
  double R_sguess;
  double psip_cguess;
  double psip_sguess;
  double psip_Rguess;

  // storing unscaled variables

  double R;
  double eta;
  double delta;
  double R_c;
  double R_s;
  double psip_c;
  double psip_s;
  double psip_R;

  // estimated upper and lower bounds for the gradient descent.
  double Rupper;
  double Rlower;
  double etaupper;
  double etalower;
  double deltaupper;
  double deltalower;
  double R_cupper;
  double R_clower;
  double R_supper;
  double R_slower;
  double psip_cupper;
  double psip_clower;
  double psip_supper;
  double psip_slower;
  double psip_Rupper;
  double psip_Rlower;
};



#endif

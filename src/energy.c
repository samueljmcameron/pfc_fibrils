#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

#define HESS(i,j) hessian[(j)+(i-1)*(3)]


// calculation of integrand arrays to be put into qromb for integration. //

void compute_rf2233b1(struct params *p,double *r,
		      double **y, double *rf_,int mpt);

void compute_integrand2(struct params *p,double *r,
			double **y,double *integrand2,int mpt);

void compute_integrand1(struct params *p,double *r,
			double **y,double *integrand1,int mpt);



// calculation of the energy (requires integrals to be calculated). //

double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt);



// calculation of first derivatives (requires integrals to be calculated). //

double derivEddelta(struct params *p,double integration2);

double derivEdeta(struct params *p,double integration1);

double derivEdR(struct params *p,double *r,double **y,
		double integration_2233b1,int mpt);




// second derivative calculations (requires integrals to be calculated). //

double ddpsidRdR(double **y,struct params *p,int mpt);

double ddEdRdR(struct params *p,double *r,double **y,
	       double integration_2233b1,int mpt);

double ddEdRdeta(struct params *p,double **y,double integration1,
		 int mpt);

double ddEdetadeta(struct params *p, double integration1);

double ddEddeltadR(struct params *p,double **y,double integration2,
		   int mpt);

double ddEddeltaddelta(struct params *p, double integration2);





void compute_rf2233b1(struct params *p,double *r,
		      double **y, double *rf_,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  rf_[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    rf_[i] = r[i]*(-(y[2][i]+0.5*sin2y/r[i])
		   +0.5*(y[2][i]+0.5*sin2y/r[i])
		   *(y[2][i]+0.5*sin2y/r[i])
		   +0.5*p->K33*siny*siny*siny
		   *siny/(r[i]*r[i])
		   +p->Lambda*p->delta*p->delta/4.0
		   //		   *(1+sin(2*eta*L)/(2*eta*L))
		   *(tmp/(cosy*cosy)-p->eta*p->eta)
		   *(tmp/(cosy*cosy)-p->eta*p->eta));
  }
  return;
}

void compute_integrand1(struct params *p,double *r,
			double **y,double *integrand1,int mpt)
{
  int i;
  double cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  integrand1[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    cosy = cos(y[1][i]);
    integrand1[i] = r[i]*(tmp/(cosy*cosy)-p->eta*p->eta);
  }
  return;
}

void compute_integrand2(struct params *p,double *r,
			double **y,double *integrand2,int mpt)
{
  int i;
  double siny, sin2y,cosy;
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);
  
  integrand2[1] = 0;
  
  for (i = 2; i <= mpt; i++) {  // compute f_fibril*r
    siny = sin(y[1][i]);
    sin2y = sin(2*y[1][i]);
    cosy = cos(y[1][i]);
    integrand2[i] = r[i]*((tmp/(cosy*cosy)-p->eta*p->eta)
			 *(tmp/(cosy*cosy)-p->eta*p->eta));
  }
  return;
}

double E_R(struct params *p,double *r,
	   double **y,double integration_2233b1,int mpt)
{

  double E;
  E = 2.0/(p->R*p->R)*integration_2233b1; // first calculate bulk energy per unit length
  // add density fluctuations term
  E = (E+p->delta*p->delta*p->omega*0.5
       *(0.75*p->delta*p->delta-1
	 //	 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
	 //	 +delta*delta*sin(4*eta*L)/(16*eta*L)
	 ));
  // add surface term tension terms
  E = E+0.5+1.0/p->R*(-(1+p->k24)*(sin(y[1][mpt])*sin(y[1][mpt]))/p->R+2.0*p->gamma_s);  
  //  E = E+2*gamma_t/L;

  return E;
}

double derivEdR(struct params *p,double *r,double **y,
		double integration_2233b1,int mpt)
{
  double ans;
  double sin2y = sin(2*y[1][mpt]);
  double siny = sin(y[1][mpt]);
  double cosy = cos(y[1][mpt]);
  double tmpcos = 4*M_PI*M_PI/(p->d0*p->d0*cosy*cosy);

  ans = -4.0/(p->R*p->R*p->R)*integration_2233b1;

  ans += 2.0/p->R*(-(y[2][mpt]+0.5*sin2y/p->R)
		   +0.5*(y[2][mpt]+0.5*sin2y/p->R)
		   *(y[2][mpt]+0.5*sin2y/p->R)
		   +0.5*p->K33*(siny*siny*siny*siny)/(p->R*p->R));

  ans += (0.5*p->Lambda*p->delta*p->delta/p->R
	  //	  *(1+sin(2*eta*L)/(2*eta*L))
	  *(tmpcos-p->eta*p->eta)
	  *(tmpcos-p->eta*p->eta));
  
  ans += 2.0/(p->R*p->R)*((1+p->k24)
			  *(siny*siny/p->R-siny*cosy*y[2][mpt])
			  -p->gamma_s);

  return ans;
}

double derivEdeta(struct params *p,double integration1)
{

  double ans=0;

  //  ans = ((0.5/(R*R)*integration2+0.5*omega*(delta*delta-1))
  //	 /eta*(cos(2*eta*L)-sin(2*eta*L)/(2*eta*L)));
  ans += (-2.0/(p->R*p->R)
	  //	  *(1+sin(2*eta*L)/(2*eta*L))
	  *p->eta*integration1);

  //  ans += (omega/(8*eta)*(cos(4*eta*L)-sin(4*eta*L)/(4*eta*L)));

  ans *= p->Lambda*p->delta*p->delta;

  return ans;
}

//double derivEdL(double Lambda,double omega,double R,double eta,
//		double delta, double gamma_t,double integration2)
//{
//  double ans;
//
//  ans = ((0.5/(R*R)*integration2
//	  +0.5*omega*(delta*delta-1))
//	 *(cos(2*eta*L)-sin(2*eta*L)/(2*eta*L)));
//
//  ans += (omega/8.0*(cos(4*eta*L)-sin(4*eta*L)/(4*eta*L)));
//  ans *= Lambda*delta*delta/L;
//  ans += -2*gamma_t/(L*L);
//
//  return ans;
//}

double derivEddelta(struct params *p,double integration2)
{
  double ans;

  ans = p->Lambda*(1.0/(p->R*p->R)
		   //	 *(1+sin(2*eta*L)/(2*eta*L))
		   *integration2);
  ans += (p->omega*p->delta*p->delta
	  *(0.75
	    //+sin(2*eta*L)/(2*eta*L)+sin(4*eta*L)/(16*eta*L)
	    ));
  ans += (p->omega*(0.75*p->delta*p->delta-1
		 //		 +sin(2*eta*L)/(2*eta*L)*(delta*delta-1)
		 //		 +delta*delta*sin(4*eta*L)/(16*eta*L)
		 ));

  ans *= p->delta;

  return ans;

}

double ddpsidRdR(double **y,struct params *p,int mpt)
{
  double ans;
  double siny = sin(y[1][mpt]);
  double sin2y = sin(2*y[1][mpt]);
  double cosy = cos(y[1][mpt]);
  double cos2y = cos(2*y[1][mpt]);
  double tmp = 4*M_PI*M_PI/(p->d0*p->d0);

  ans = (1.0/p->R*(1-y[2][mpt]
		   +p->K33*siny*siny*sin2y/p->R
		   -cos2y*(1-0.5*sin2y/p->R)
		   +tmp*p->Lambda*p->delta*p->delta
		   *(tmp/(cosy*cosy)-p->eta*p->eta)
		   *siny*p->R/(cosy*cosy*cosy)));
  return ans;
}

double ddEdRdR(struct params *p,double *r,double **y,
	       double integration_2233b1,int mpt)
{
  double ans;

  double sin2y = sin(2*y[1][mpt]);
  double siny = sin(y[1][mpt]);
  double cosy = cos(y[1][mpt]);
  double cos2y = cos(2*y[1][mpt]);
  double tmpcos = 4*M_PI*M_PI/(p->d0*p->d0*cosy*cosy);

  ans = 12.0/(p->R*p->R*p->R*p->R)*integration_2233b1;

  ans += -8.0/(p->R*p->R)*(-(y[2][mpt]+0.5*sin2y/p->R)
			   +0.5*(y[2][mpt]+0.5*sin2y/p->R)
			   *(y[2][mpt]+0.5*sin2y/p->R)
			   +0.5*p->K33*(siny*siny*siny*siny)/(p->R*p->R));
  ans += 2.0/p->R*((ddpsidRdR(y,p,mpt)-0.5*sin2y/(p->R*p->R)
		    +cos2y*y[2][mpt]/p->R)*(y[2][mpt]+0.5*sin2y/p->R-1)
		   +2*p->K33*siny*siny*siny*cosy*y[2][mpt]/(p->R*p->R));

  ans += (-0.5*p->Lambda*p->delta*p->delta/(p->R*p->R)
	  //	  *(1+sin(2*eta*L)/(2*eta*L))
	  *(tmpcos-p->eta*p->eta)
	  *(tmpcos-p->eta*p->eta));

  ans += (p->Lambda*p->delta*p->delta/p->R*(tmpcos-p->eta*p->eta)
	  *2*tmpcos*siny/cosy*y[2][mpt]);

  ans += -6.0*(1+p->k24)/(p->R*p->R*p->R*p->R)*siny*siny;

  ans += 4.0*(1+p->k24)/(p->R*p->R*p->R)*siny*cosy*y[2][mpt];

  ans += 2.0*(1+p->k24)/(p->R*p->R*p->R)*sin2y*y[2][mpt];

  ans += -2.0*(1+p->k24)/(p->R*p->R)*cos2y*y[2][mpt]*y[2][mpt];

  ans += -1.0*(1+p->k24)/(p->R*p->R)*sin2y*ddpsidRdR(y,p,mpt);

  ans += 4.0/(p->R*p->R*p->R)*p->gamma_s;

  return ans;

}

double ddEdRdeta(struct params *p,double **y,double integration1,
		 int mpt)
{
  double ans;
  double cosy = cos(y[1][mpt]);

  ans = 4/(p->R*p->R*p->R)*integration1;

  ans += -2.0/p->R*(4*M_PI*M_PI/(p->d0*p->d0*cosy*cosy)-p->eta*p->eta);

  ans *= p->Lambda*p->delta*p->delta*p->eta;

  return ans;
}

double ddEdetadeta(struct params *p, double integration1)
{

  double ans;

  ans = -2.0/(p->R*p->R)*integration1+2.0*p->eta;

  ans *= p->Lambda*p->delta*p->delta;

  return ans;
}

double ddEddeltadR(struct params *p,double **y,double integration2,
		   int mpt)
{

  double ans;
  double cosy = cos(y[1][mpt]);
  double tmpcos = 4*M_PI*M_PI/(p->d0*p->d0*cosy*cosy);

  ans = -2.0/(p->R*p->R*p->R)*integration2;
  
  ans += 1.0/p->R*(tmpcos-p->eta*p->eta)*(tmpcos-p->eta*p->eta);

  ans *= p->Lambda*p->delta;

  return ans;
}

double ddEddeltadeta(struct params *p, double integration1)
{

  double ans;

  ans = -4.0/(p->R*p->R)*integration1;
  
  ans *= p->Lambda*p->delta*p->eta;

  return ans;
}

double ddEddeltaddelta(struct params *p, double integration2)
{

  double ans;

  ans = p->Lambda/(p->R*p->R)*integration2;

  ans += p->omega*(4.5*p->delta*p->delta-1);

  return ans;
}

bool energy_properties(double *E, double *dEdx,struct params *p,
		       double *r,double **y,double *rf_,
		       double *integrand1,double *integrand2,int mpt)
{
  bool failure = true;
  double integration_2233b1,integration1,integration2;
  double tol0 = 1e-14;
  double tol2233b1,tol1,tol2;

  // tolerances ("tol...") are to determine how large the energy or 
  // derivative terms are without the integrals (setting them to 0).
  // The reason for this is if the integral calculation performed
  // by qromb has not converged (usually because the integral error
  // is so small that round-off error becomes an issue), if the
  // size of the error term is so small relative to the rest of the
  // function that it doesn't change the functions value (up to
  // some tolerance tolblabla), then the lack of convergence can just
  // be ignored.
  
  tol2233b1 = fabs(E_R(p,r,y,0,mpt)*p->R*p->R/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,r,y,rf_,mpt);
  integration_2233b1 = qromb(r,rf_,mpt,tol2233b1,&failure);

  if (failure) {
    printf("tol2233b1=%e\n",tol2233b1);
    //    printf("failure occurred at integration_2233b1.\n");
    return false;
  }

  if (dEdx[1] != 0 || dEdx[2] != 0 || dEdx[3] != 0) {
        
    if (fabs(p->delta)<=tol0) {
      dEdx[2] = 0;
      
      dEdx[3] = 0;


    }
    else if (fabs(p->eta) <= tol0) {
      dEdx[2] = 0;

      tol2 = fabs(derivEddelta(p,0)*p->R*p->R/(p->Lambda*p->delta)*tol0);
      tol2 = tol2 > tol0 ? tol2 : tol0;
      compute_integrand2(p,r,y,integrand2,mpt);
      integration2 = qromb(r,integrand2,mpt,tol2,&failure);
      
      if (failure) {
	printf("tol2=%e\n",tol2);
	//    printf("failure occurred at integration2.\n");
	return false;
      }


      dEdx[3] = derivEddelta(p,integration2);


    }
    else {
      

      tol1 = tol0;
      compute_integrand1(p,r,y,integrand1,mpt);
      integration1 = qromb(r,integrand1,mpt,tol1,&failure);
      
      if (failure) {
	printf("tol1=%e\n",tol1);
	//    printf("failure occurred at integration1.\n");
	return false;
      }
      
      tol2 = fabs(derivEddelta(p,0)*p->R*p->R/(p->Lambda*p->delta)*tol0);
      tol2 = tol2 > tol0 ? tol2 : tol0;
      compute_integrand2(p,r,y,integrand2,mpt);
      integration2 = qromb(r,integrand2,mpt,tol2,&failure);
      
      if (failure) {
	printf("tol2=%e\n",tol2);
	//    printf("failure occurred at integration2.\n");
	return false;
      }
      
      dEdx[2] = derivEdeta(p,integration1);
      
      dEdx[3] = derivEddelta(p,integration2);


    }
    
    dEdx[1] = derivEdR(p,r,y,integration_2233b1,mpt);

    // the second derivatives with respect to R, and with respect to delta,
    // are both non-zero regardless of parameter values.

  }
  *E = E_R(p,r,y,integration_2233b1,mpt);

  return true;
}


bool energy_prop_with_hessian(double *E, double *dEdx,struct params *p,
			      double *r,double **y,double *rf_,
			      double *integrand1,double *integrand2,int mpt,
			      double *hessian)
{
  bool failure = true;
  double integration_2233b1,integration1,integration2;
  double tol0 = 1e-14;
  double tol2233b1,tol1,tol2;

  // tolerances ("tol...") are to determine how large the energy or 
  // derivative terms are without the integrals (setting them to 0).
  // The reason for this is if the integral calculation performed
  // by qromb has not converged (usually because the integral error
  // is so small that round-off error becomes an issue), if the
  // size of the error term is so small relative to the rest of the
  // function that it doesn't change the functions value (up to
  // some tolerance tolblabla), then the lack of convergence can just
  // be ignored.
  
  tol2233b1 = fabs(E_R(p,r,y,0,mpt)*p->R*p->R/2.0*tol0);
  tol2233b1 = tol2233b1 > tol0 ? tol2233b1 : tol0;
  compute_rf2233b1(p,r,y,rf_,mpt);
  integration_2233b1 = qromb(r,rf_,mpt,tol2233b1,&failure);

  if (failure) {
    printf("tol2233b1=%e\n",tol2233b1);
    //    printf("failure occurred at integration_2233b1.\n");
    return false;
  }

  if (dEdx[1] != 0 || dEdx[2] != 0 || dEdx[3] != 0) {
    
    tol2 = fabs(derivEddelta(p,0)*p->R*p->R/(p->Lambda*p->delta)*tol0);
    tol2 = tol2 > tol0 ? tol2 : tol0;
    compute_integrand2(p,r,y,integrand2,mpt);
    integration2 = qromb(r,integrand2,mpt,tol2,&failure);
    
    if (failure) {
      printf("tol2=%e\n",tol2);
      //    printf("failure occurred at integration2.\n");
      return false;
    }
    
    if (fabs(p->delta)<=tol0) {
      dEdx[2] = 0;
      
      dEdx[3] = 0;

      HESS(1,2) = HESS(2,1) = 0;

      HESS(1,3) = HESS(3,1) = 0;

      HESS(2,3) = HESS(3,2) = 0;

      HESS(2,2) = 0;

    }
    else if (fabs(p->eta) <= tol0) {
      dEdx[2] = 0;
      
      HESS(1,2) = HESS(2,1) = 0;

      HESS(2,3) = HESS(3,2) = 0;

      tol1 = tol0;
      compute_integrand1(p,r,y,integrand1,mpt);
      integration1 = qromb(r,integrand1,mpt,tol1,&failure);
      
      if (failure) {
	printf("tol1=%e\n",tol1);
	//    printf("failure occurred at integration1.\n");
	return false;
      }

      dEdx[3] = derivEddelta(p,integration2);

      HESS(1,3) = HESS(3,1) = ddEddeltadR(p,y,integration2,mpt);

      HESS(2,2) = ddEdetadeta(p,integration1);

    }
    else {
      

      tol1 = tol0;
      compute_integrand1(p,r,y,integrand1,mpt);
      integration1 = qromb(r,integrand1,mpt,tol1,&failure);
      
      if (failure) {
	printf("tol1=%e\n",tol1);
	//    printf("failure occurred at integration1.\n");
	return false;
      }
      
      tol2 = fabs(derivEddelta(p,0)*p->R*p->R/(p->Lambda*p->delta)*tol0);
      tol2 = tol2 > tol0 ? tol2 : tol0;
      compute_integrand2(p,r,y,integrand2,mpt);
      integration2 = qromb(r,integrand2,mpt,tol2,&failure);
      
      if (failure) {
	printf("tol2=%e\n",tol2);
	//    printf("failure occurred at integration2.\n");
	return false;
      }
      
      dEdx[2] = derivEdeta(p,integration1);
      
      dEdx[3] = derivEddelta(p,integration2);

      HESS(1,2) = HESS(2,1) = ddEdRdeta(p,y,integration1,mpt);

      HESS(1,3) = HESS(3,1) = ddEddeltadR(p,y,integration2,mpt);

      HESS(2,3) = HESS(3,2) = ddEddeltadeta(p,integration1);

      HESS(2,2) = ddEdetadeta(p,integration1);

    }
    
    dEdx[1] = derivEdR(p,r,y,integration_2233b1,mpt);

    // the second derivatives with respect to R, and with respect to delta,
    // are both non-zero regardless of parameter values.

    HESS(1,1) = ddEdRdR(p,r,y,integration_2233b1,mpt);

    HESS(3,3) = ddEddeltaddelta(p,integration2);

  }
  *E = E_R(p,r,y,integration_2233b1,mpt);

  return true;
}


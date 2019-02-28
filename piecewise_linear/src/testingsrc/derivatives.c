#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include "edited_gsl_src/gsl_multimin.h"
#include "funcs/nrutil.h"
#include "headerfile.h"

double Efunc(struct params *p);
double dEdR(struct params *p);
double dEdeta(struct params *p);
double dEddelta(struct params *p);
double dEdR_c(struct params *p);
double dEdR_s(struct params *p);
double dEdpsip_c(struct params *p);
double dEdpsip_s(struct params *p);
double dEdpsip_R(struct params *p);


int main(int argc, char **argv)
{

  void setparams(struct params *p,char **args);
  
  void test_dEdR(struct params p0,struct params p1,struct params p2);

  void test_dEdeta(struct params p0,struct params p1,struct params p2);

  void test_dEddelta(struct params p0,struct params p1,struct params p2);

  void test_dEdR_c(struct params p0,struct params p1,struct params p2);
  
  void test_dEdR_s(struct params p0,struct params p1,struct params p2);
  
  void test_dEdpsip_c(struct params p0,struct params p1,struct params p2);
  
  void test_dEdpsip_s(struct params p0,struct params p1,struct params p2);
  
  void test_dEdpsip_R(struct params p0,struct params p1,struct params p2);

  struct params p0,p1,p2;

  setparams(&p0,argv);

  setparams(&p1,argv);
  setparams(&p2,argv);

  printf("E1 = %lf\n",Efunc(&p0));
  
  test_dEdR(p0,p1,p2);
  test_dEdeta(p0,p1,p2);
  test_dEddelta(p0,p1,p2);
  test_dEdR_c(p0,p1,p2);
  test_dEdR_s(p0,p1,p2);

  test_dEdpsip_c(p0,p1,p2);

  test_dEdpsip_s(p0,p1,p2);

  test_dEdpsip_R(p0,p1,p2);


  return 0;

}

void setparams(struct params *p,char **args)
{
  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->R);
  sscanf(args[8],"%lf",&p->eta);
  sscanf(args[9],"%lf",&p->delta);
  sscanf(args[10],"%lf",&p->R_c);
  sscanf(args[11],"%lf",&p->R_s);
  sscanf(args[12],"%lf",&p->psip_c);
  sscanf(args[13],"%lfn",&p->psip_s);
  sscanf(args[14],"%lf",&p->psip_R);
  p->psip_s = 0.0005;
  
  return;
}


void test_dEdR(struct params p0,struct params p1,struct params p2)
{
  
  double dEdR_i;
  double dEdR_ii;
  double dEdR_iii;
  double dEdR_iv;
  double dEdR_v;
  double dEdR_vi;
  double dEdR_vii;
  
  double exact_dEdR;
  
  double dR_i = 1e-2;
  double dR_ii = 1e-3;
  double dR_iii = 1e-4;
  double dR_iv = 1e-5;
  double dR_v = 1e-6;
  double dR_vi = 1e-7;
  double dR_vii = 1e-8;


  p1.R = p0.R-dR_i;
  p2.R = p0.R+dR_i;
  dEdR_i = (Efunc(&p2)-Efunc(&p1))/(2*dR_i);

  p1.R = p0.R-dR_ii;
  p2.R = p0.R+dR_ii;
  dEdR_ii = (Efunc(&p2)-Efunc(&p1))/(2*dR_ii);

  p1.R = p0.R-dR_iii;
  p2.R = p0.R+dR_iii;
  dEdR_iii = (Efunc(&p2)-Efunc(&p1))/(2*dR_iii);

  p1.R = p0.R-dR_iv;
  p2.R = p0.R+dR_iv;
  dEdR_iv = (Efunc(&p2)-Efunc(&p1))/(2*dR_iv);
  
  p1.R = p0.R-dR_v;
  p2.R = p0.R+dR_v;
  dEdR_v = (Efunc(&p2)-Efunc(&p1))/(2*dR_v);

  p1.R = p0.R-dR_vi;
  p2.R = p0.R+dR_vi;
  dEdR_vi = (Efunc(&p2)-Efunc(&p1))/(2*dR_vi);

  p1.R = p0.R-dR_vii;
  p2.R = p0.R+dR_vii;
  dEdR_vii = (Efunc(&p2)-Efunc(&p1))/(2*dR_vii);
  
  exact_dEdR = dEdR(&p0);

  printf("dEdR = %.6e, dR = %e\n",dEdR_i,dR_i);
  printf("dEdR = %.6e, dR = %e\n",dEdR_ii,dR_ii);
  printf("dEdR = %.6e, dR = %e\n",dEdR_iii,dR_iii);
  printf("dEdR = %.6e, dR = %e\n",dEdR_iv,dR_iv);
  printf("dEdR = %.6e, dR = %e\n",dEdR_v,dR_v);
  printf("dEdR = %.6e, dR = %e\n",dEdR_vi,dR_vi);
  printf("dEdR = %.6e, dR = %e\n",dEdR_vii,dR_vii);
  
  printf("dEdR = %.6e, dR = 0\n",exact_dEdR);


  return;

}


void test_dEdeta(struct params p0,struct params p1,struct params p2)
{
  
  double dEdeta_i;
  double dEdeta_ii;
  double dEdeta_iii;
  double dEdeta_iv;
  double dEdeta_v;
  double dEdeta_vi;
  double dEdeta_vii;
  
  double exact_dEdeta;
  
  double deta_i = 1e-2;
  double deta_ii = 1e-3;
  double deta_iii = 1e-4;
  double deta_iv = 1e-5;
  double deta_v = 1e-6;
  double deta_vi = 1e-7;
  double deta_vii = 1e-8;


  p1.eta = p0.eta-deta_i;
  p2.eta = p0.eta+deta_i;
  dEdeta_i = (Efunc(&p2)-Efunc(&p1))/(2*deta_i);

  p1.eta = p0.eta-deta_ii;
  p2.eta = p0.eta+deta_ii;
  dEdeta_ii = (Efunc(&p2)-Efunc(&p1))/(2*deta_ii);

  p1.eta = p0.eta-deta_iii;
  p2.eta = p0.eta+deta_iii;
  dEdeta_iii = (Efunc(&p2)-Efunc(&p1))/(2*deta_iii);

  p1.eta = p0.eta-deta_iv;
  p2.eta = p0.eta+deta_iv;
  dEdeta_iv = (Efunc(&p2)-Efunc(&p1))/(2*deta_iv);
  
  p1.eta = p0.eta-deta_v;
  p2.eta = p0.eta+deta_v;
  dEdeta_v = (Efunc(&p2)-Efunc(&p1))/(2*deta_v);

  p1.eta = p0.eta-deta_vi;
  p2.eta = p0.eta+deta_vi;
  dEdeta_vi = (Efunc(&p2)-Efunc(&p1))/(2*deta_vi);

  p1.eta = p0.eta-deta_vii;
  p2.eta = p0.eta+deta_vii;
  dEdeta_vii = (Efunc(&p2)-Efunc(&p1))/(2*deta_vii);
  
  exact_dEdeta = dEdeta(&p0);

  printf("dEdeta = %.6e, deta = %e\n",dEdeta_i,deta_i);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_ii,deta_ii);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_iii,deta_iii);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_iv,deta_iv);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_v,deta_v);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_vi,deta_vi);
  printf("dEdeta = %.6e, deta = %e\n",dEdeta_vii,deta_vii);
  
  printf("dEdeta = %.6e, deta = 0\n",exact_dEdeta);


  return;

}


void test_dEddelta(struct params p0,struct params p1,struct params p2)
{
  
  double dEddelta_i;
  double dEddelta_ii;
  double dEddelta_iii;
  double dEddelta_iv;
  double dEddelta_v;
  double dEddelta_vi;
  double dEddelta_vii;
  
  double exact_dEddelta;
  
  double ddelta_i = 1e-2;
  double ddelta_ii = 1e-3;
  double ddelta_iii = 1e-4;
  double ddelta_iv = 1e-5;
  double ddelta_v = 1e-6;
  double ddelta_vi = 1e-7;
  double ddelta_vii = 1e-8;


  p1.delta = p0.delta-ddelta_i;
  p2.delta = p0.delta+ddelta_i;
  dEddelta_i = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_i);

  p1.delta = p0.delta-ddelta_ii;
  p2.delta = p0.delta+ddelta_ii;
  dEddelta_ii = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_ii);

  p1.delta = p0.delta-ddelta_iii;
  p2.delta = p0.delta+ddelta_iii;
  dEddelta_iii = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_iii);

  p1.delta = p0.delta-ddelta_iv;
  p2.delta = p0.delta+ddelta_iv;
  dEddelta_iv = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_iv);
  
  p1.delta = p0.delta-ddelta_v;
  p2.delta = p0.delta+ddelta_v;
  dEddelta_v = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_v);

  p1.delta = p0.delta-ddelta_vi;
  p2.delta = p0.delta+ddelta_vi;
  dEddelta_vi = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_vi);

  p1.delta = p0.delta-ddelta_vii;
  p2.delta = p0.delta+ddelta_vii;
  dEddelta_vii = (Efunc(&p2)-Efunc(&p1))/(2*ddelta_vii);
  
  exact_dEddelta = dEddelta(&p0);

  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_i,ddelta_i);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_ii,ddelta_ii);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_iii,ddelta_iii);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_iv,ddelta_iv);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_v,ddelta_v);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_vi,ddelta_vi);
  printf("dEddelta = %.6e, ddelta = %e\n",dEddelta_vii,ddelta_vii);
  
  printf("dEddelta = %.6e, ddelta = 0\n",exact_dEddelta);


  return;

}

void test_dEdR_c(struct params p0,struct params p1,struct params p2)
{
  
  double dEdR_c_i;
  double dEdR_c_ii;
  double dEdR_c_iii;
  double dEdR_c_iv;
  double dEdR_c_v;
  double dEdR_c_vi;
  double dEdR_c_vii;
  
  double exact_dEdR_c;
  
  double dR_c_i = 1e-2;
  double dR_c_ii = 1e-3;
  double dR_c_iii = 1e-4;
  double dR_c_iv = 1e-5;
  double dR_c_v = 1e-6;
  double dR_c_vi = 1e-7;
  double dR_c_vii = 1e-8;


  p1.R_c = p0.R_c-dR_c_i;
  p2.R_c = p0.R_c+dR_c_i;
  dEdR_c_i = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_i);

  p1.R_c = p0.R_c-dR_c_ii;
  p2.R_c = p0.R_c+dR_c_ii;
  dEdR_c_ii = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_ii);

  p1.R_c = p0.R_c-dR_c_iii;
  p2.R_c = p0.R_c+dR_c_iii;
  dEdR_c_iii = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_iii);

  p1.R_c = p0.R_c-dR_c_iv;
  p2.R_c = p0.R_c+dR_c_iv;
  dEdR_c_iv = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_iv);
  
  p1.R_c = p0.R_c-dR_c_v;
  p2.R_c = p0.R_c+dR_c_v;
  dEdR_c_v = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_v);

  p1.R_c = p0.R_c-dR_c_vi;
  p2.R_c = p0.R_c+dR_c_vi;
  dEdR_c_vi = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_vi);

  p1.R_c = p0.R_c-dR_c_vii;
  p2.R_c = p0.R_c+dR_c_vii;
  dEdR_c_vii = (Efunc(&p2)-Efunc(&p1))/(2*dR_c_vii);
  
  exact_dEdR_c = dEdR_c(&p0);

  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_i,dR_c_i);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_ii,dR_c_ii);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_iii,dR_c_iii);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_iv,dR_c_iv);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_v,dR_c_v);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_vi,dR_c_vi);
  printf("dEdR_c = %.6e, dR_c = %e\n",dEdR_c_vii,dR_c_vii);
  
  printf("dEdR_c = %.6e, dR_c = 0\n",exact_dEdR_c);


  return;

}


void test_dEdR_s(struct params p0,struct params p1,struct params p2)
{
  
  double dEdR_s_i;
  double dEdR_s_ii;
  double dEdR_s_iii;
  double dEdR_s_iv;
  double dEdR_s_v;
  double dEdR_s_vi;
  double dEdR_s_vii;
  
  double exact_dEdR_s;
  
  double dR_s_i = 1e-2;
  double dR_s_ii = 1e-3;
  double dR_s_iii = 1e-4;
  double dR_s_iv = 1e-5;
  double dR_s_v = 1e-6;
  double dR_s_vi = 1e-7;
  double dR_s_vii = 1e-8;


  p1.R_s = p0.R_s-dR_s_i;
  p2.R_s = p0.R_s+dR_s_i;
  dEdR_s_i = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_i);

  p1.R_s = p0.R_s-dR_s_ii;
  p2.R_s = p0.R_s+dR_s_ii;
  dEdR_s_ii = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_ii);

  p1.R_s = p0.R_s-dR_s_iii;
  p2.R_s = p0.R_s+dR_s_iii;
  dEdR_s_iii = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_iii);

  p1.R_s = p0.R_s-dR_s_iv;
  p2.R_s = p0.R_s+dR_s_iv;
  dEdR_s_iv = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_iv);
  
  p1.R_s = p0.R_s-dR_s_v;
  p2.R_s = p0.R_s+dR_s_v;
  dEdR_s_v = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_v);

  p1.R_s = p0.R_s-dR_s_vi;
  p2.R_s = p0.R_s+dR_s_vi;
  dEdR_s_vi = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_vi);

  p1.R_s = p0.R_s-dR_s_vii;
  p2.R_s = p0.R_s+dR_s_vii;
  dEdR_s_vii = (Efunc(&p2)-Efunc(&p1))/(2*dR_s_vii);
  
  exact_dEdR_s = dEdR_s(&p0);

  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_i,dR_s_i);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_ii,dR_s_ii);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_iii,dR_s_iii);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_iv,dR_s_iv);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_v,dR_s_v);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_vi,dR_s_vi);
  printf("dEdR_s = %.6e, dR_s = %e\n",dEdR_s_vii,dR_s_vii);
  
  printf("dEdR_s = %.6e, dR_s = 0\n",exact_dEdR_s);


  return;

}

void test_dEdpsip_c(struct params p0,struct params p1,struct params p2)
{
  
  double dEdpsip_c_i;
  double dEdpsip_c_ii;
  double dEdpsip_c_iii;
  double dEdpsip_c_iv;
  double dEdpsip_c_v;
  double dEdpsip_c_vi;
  double dEdpsip_c_vii;
  
  double exact_dEdpsip_c;
  
  double dpsip_c_i = 1e-2;
  double dpsip_c_ii = 1e-3;
  double dpsip_c_iii = 1e-4;
  double dpsip_c_iv = 1e-5;
  double dpsip_c_v = 1e-6;
  double dpsip_c_vi = 1e-7;
  double dpsip_c_vii = 1e-8;


  p1.psip_c = p0.psip_c-dpsip_c_i;
  p2.psip_c = p0.psip_c+dpsip_c_i;
  dEdpsip_c_i = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_i);

  p1.psip_c = p0.psip_c-dpsip_c_ii;
  p2.psip_c = p0.psip_c+dpsip_c_ii;
  dEdpsip_c_ii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_ii);

  p1.psip_c = p0.psip_c-dpsip_c_iii;
  p2.psip_c = p0.psip_c+dpsip_c_iii;
  dEdpsip_c_iii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_iii);

  p1.psip_c = p0.psip_c-dpsip_c_iv;
  p2.psip_c = p0.psip_c+dpsip_c_iv;
  dEdpsip_c_iv = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_iv);
  
  p1.psip_c = p0.psip_c-dpsip_c_v;
  p2.psip_c = p0.psip_c+dpsip_c_v;
  dEdpsip_c_v = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_v);

  p1.psip_c = p0.psip_c-dpsip_c_vi;
  p2.psip_c = p0.psip_c+dpsip_c_vi;
  dEdpsip_c_vi = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_vi);

  p1.psip_c = p0.psip_c-dpsip_c_vii;
  p2.psip_c = p0.psip_c+dpsip_c_vii;
  dEdpsip_c_vii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_c_vii);
  
  exact_dEdpsip_c = dEdpsip_c(&p0);

  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_i,dpsip_c_i);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_ii,dpsip_c_ii);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_iii,dpsip_c_iii);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_iv,dpsip_c_iv);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_v,dpsip_c_v);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_vi,dpsip_c_vi);
  printf("dEdpsip_c = %.6e, dpsip_c = %e\n",dEdpsip_c_vii,dpsip_c_vii);
  
  printf("dEdpsip_c = %.6e, dpsip_c = 0\n",exact_dEdpsip_c);


  return;

}

void test_dEdpsip_s(struct params p0,struct params p1,struct params p2)
{
  
  double dEdpsip_s_i;
  double dEdpsip_s_ii;
  double dEdpsip_s_iii;
  double dEdpsip_s_iv;
  double dEdpsip_s_v;
  double dEdpsip_s_vi;
  double dEdpsip_s_vii;
  
  double exact_dEdpsip_s;
  
  double dpsip_s_i = 1e-2;
  double dpsip_s_ii = 1e-3;
  double dpsip_s_iii = 1e-4;
  double dpsip_s_iv = 1e-5;
  double dpsip_s_v = 1e-6;
  double dpsip_s_vi = 1e-7;
  double dpsip_s_vii = 1e-8;


  p1.psip_s = p0.psip_s-dpsip_s_i;
  p2.psip_s = p0.psip_s+dpsip_s_i;
  dEdpsip_s_i = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_i);

  p1.psip_s = p0.psip_s-dpsip_s_ii;
  p2.psip_s = p0.psip_s+dpsip_s_ii;
  dEdpsip_s_ii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_ii);

  p1.psip_s = p0.psip_s-dpsip_s_iii;
  p2.psip_s = p0.psip_s+dpsip_s_iii;
  dEdpsip_s_iii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_iii);

  p1.psip_s = p0.psip_s-dpsip_s_iv;
  p2.psip_s = p0.psip_s+dpsip_s_iv;
  dEdpsip_s_iv = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_iv);
  
  p1.psip_s = p0.psip_s-dpsip_s_v;
  p2.psip_s = p0.psip_s+dpsip_s_v;
  dEdpsip_s_v = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_v);

  p1.psip_s = p0.psip_s-dpsip_s_vi;
  p2.psip_s = p0.psip_s+dpsip_s_vi;
  dEdpsip_s_vi = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_vi);

  p1.psip_s = p0.psip_s-dpsip_s_vii;
  p2.psip_s = p0.psip_s+dpsip_s_vii;
  dEdpsip_s_vii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_s_vii);
  
  exact_dEdpsip_s = dEdpsip_s(&p0);

  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_i,dpsip_s_i);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_ii,dpsip_s_ii);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_iii,dpsip_s_iii);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_iv,dpsip_s_iv);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_v,dpsip_s_v);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_vi,dpsip_s_vi);
  printf("dEdpsip_s = %.6e, dpsip_s = %e\n",dEdpsip_s_vii,dpsip_s_vii);
  
  printf("dEdpsip_s = %.6e, dpsip_s = 0\n",exact_dEdpsip_s);


  return;

}

void test_dEdpsip_R(struct params p0,struct params p1,struct params p2)
{
  
  double dEdpsip_R_i;
  double dEdpsip_R_ii;
  double dEdpsip_R_iii;
  double dEdpsip_R_iv;
  double dEdpsip_R_v;
  double dEdpsip_R_vi;
  double dEdpsip_R_vii;
  
  double exact_dEdpsip_R;
  
  double dpsip_R_i = 1e-2;
  double dpsip_R_ii = 1e-3;
  double dpsip_R_iii = 1e-4;
  double dpsip_R_iv = 1e-5;
  double dpsip_R_v = 1e-6;
  double dpsip_R_vi = 1e-7;
  double dpsip_R_vii = 1e-8;


  p1.psip_R = p0.psip_R-dpsip_R_i;
  p2.psip_R = p0.psip_R+dpsip_R_i;
  dEdpsip_R_i = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_i);

  p1.psip_R = p0.psip_R-dpsip_R_ii;
  p2.psip_R = p0.psip_R+dpsip_R_ii;
  dEdpsip_R_ii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_ii);

  p1.psip_R = p0.psip_R-dpsip_R_iii;
  p2.psip_R = p0.psip_R+dpsip_R_iii;
  dEdpsip_R_iii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_iii);

  p1.psip_R = p0.psip_R-dpsip_R_iv;
  p2.psip_R = p0.psip_R+dpsip_R_iv;
  dEdpsip_R_iv = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_iv);
  
  p1.psip_R = p0.psip_R-dpsip_R_v;
  p2.psip_R = p0.psip_R+dpsip_R_v;
  dEdpsip_R_v = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_v);

  p1.psip_R = p0.psip_R-dpsip_R_vi;
  p2.psip_R = p0.psip_R+dpsip_R_vi;
  dEdpsip_R_vi = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_vi);

  p1.psip_R = p0.psip_R-dpsip_R_vii;
  p2.psip_R = p0.psip_R+dpsip_R_vii;
  dEdpsip_R_vii = (Efunc(&p2)-Efunc(&p1))/(2*dpsip_R_vii);
  
  exact_dEdpsip_R = dEdpsip_R(&p0);

  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_i,dpsip_R_i);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_ii,dpsip_R_ii);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_iii,dpsip_R_iii);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_iv,dpsip_R_iv);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_v,dpsip_R_v);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_vi,dpsip_R_vi);
  printf("dEdpsip_R = %.6e, dpsip_R = %e\n",dEdpsip_R_vii,dpsip_R_vii);
  
  printf("dEdpsip_R = %.6e, dpsip_R = 0\n",exact_dEdpsip_R);


  return;

}

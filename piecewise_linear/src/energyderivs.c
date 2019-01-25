#include <stdio.h>
#include <stdlib.h>
#include "headerfile.h"


double dudx_1(double x_1,double zeta);
double dudx_2(double x_2,double zeta);
double dudzeta(double x_1,double x_2,double zeta);
double dvdx_1(double x_1,double xi,double zeta);
double dvdx_2(double x_2,double xi,double zeta);
double dvdxi(double x_1,double x_2,double xi, double zeta);
double dvdzeta(double x_1,double x_2,double xi,double zeta);
double df_1dx_1(double x_1,double xi,double zeta);
double df_1dx_2(double x_2,double xi,double zeta);
double df_1dxi(double x_1,double x_2,double xi,double zeta);
double df_2dxi(double x_1,double x_2,double xi,double zeta);
double df_1dzeta(double x_1,double x_2,double xi,double zeta);
double df_2dzeta(double x_1,double x_2,double xi,double zeta);
double dg_1dx_1(double x_1,double xi,double zeta);
double dg_1dx_2(double x_2,double xi,double zeta);
double dg_1dxi(double x_1,double x_2,double xi,double zeta);
double dg_2dxi(double x_1,double x_2,double xi,double zeta);
double dg_1dzeta(double x_1,double x_2,double xi,double zeta);
double dg_2dzeta(double x_1,double x_2,double xi,double zeta);

double dEdR(double *x, struct params *p)
{



}

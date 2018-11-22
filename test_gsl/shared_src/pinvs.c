#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "nrutil.h"
#include "headerfile.h"

bool pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, double ***c,
	   double **s)
/*Diagonalize the square subsection of the s matrix, and store the recursion coefficients in c; used internally by solvde.*/
{
  int js1,jpiv,jp,je2,jcoff,j,irow,ipiv,id,icoff,i,*indxr;
  double pivinv,piv,dum,big,*pscl;
  indxr=ivector(ie1,ie2);
  pscl=vector(ie1,ie2);
  je2=je1+ie2-ie1;
  js1=je2+1;
  for (i=ie1;i<=ie2;i++) {// Implicit pivoting, as in §2.1.
    big=0.0;
    for (j=je1;j<=je2;j++)
      if (fabs(s[i][j]) > big) big=fabs(s[i][j]);
    if (big == 0.0) {
      printf("Singular matrix - row all 0, in pinvs\n");
      return false;
    }
    pscl[i]=1.0/big;
    indxr[i]=0;
  }
  for (id=ie1;id<=ie2;id++) {
    piv=0.0;
    for (i=ie1;i<=ie2;i++) { //Find pivot element.
      if (indxr[i] == 0) {
	big=0.0;
	for (j=je1;j<=je2;j++) {
	  if (fabs(s[i][j]) > big) {
	    jp=j;
	    big=fabs(s[i][j]);
	  }
	}
	if (big*pscl[i] > piv) {
	  ipiv=i;
	  jpiv=jp;
	  piv=big*pscl[i];
	}
      }
    }
    if (s[ipiv][jpiv] == 0.0) {
      printf("Singular matrix in routine pinvs\n");
      return false;
    }
    indxr[ipiv]=jpiv; //In place reduction. Save column ordering.
    pivinv=1.0/s[ipiv][jpiv];
    for (j=je1;j<=jsf;j++) s[ipiv][j] *= pivinv; //Normalize pivot row.
    s[ipiv][jpiv]=1.0;
    for (i=ie1;i<=ie2;i++) { //Reduce nonpivot elements in column.
      if (indxr[i] != jpiv) {
	if (s[i][jpiv]) {
	  dum=s[i][jpiv];
	  for (j=je1;j<=jsf;j++)
	    s[i][j] -= dum*s[ipiv][j];
	  s[i][jpiv]=0.0;
	}
      }
    }
  }
  jcoff=jc1-js1; //Sort and store unreduced coefficients.
  icoff=ie1-je1;
  for (i=ie1;i<=ie2;i++) {
    irow=indxr[i]+icoff;
    for (j=js1;j<=jsf;j++) c[irow][j+jcoff][k]=s[i][j];
  }
  free_vector(pscl,ie1,ie2);
  free_ivector(indxr,ie1,ie2);
  return true;
}

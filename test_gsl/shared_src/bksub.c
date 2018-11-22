#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"


void bksub(int ne, int nb, int jf, int k1, int k2, double ***c)
//Backsubstitution, used internally by solvde.
{
  int nbf,im,kp,k,j,i;
  double xx;
  nbf=ne-nb;
  im=1;
  for (k=k2;k>=k1;k--) { //Use recurrence relations to eliminate remaining de
    if (k == k1) im=nbf+1;   // pendences.
    kp=k+1;
    for (j=1;j<=nbf;j++) {
      xx=c[j][jf][kp];
      for (i=im;i<=ne;i++)
	c[i][jf][k] -= c[i][j][k]*xx;
    }
  }
  for (k=k1;k<=k2;k++) { //Reorder corrections to be in column 1.
    kp=k+1;
    for (i=1;i<=nb;i++) c[i][1][k]=c[i+nbf][jf][k];
    for (i=1;i<=nbf;i++) c[i+nb][1][k]=c[i][jf][kp];
  }
}

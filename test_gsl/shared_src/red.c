#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "headerfile.h"

void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	 int ic1, int jc1, int jcf, int kc, double ***c, double **s)
/*Reduce columns jz1-jz2 of the s matrix, using previous results as stored in the c matrix. Only columns jm1-jm2,jmf are affected by the prior results. red is used internally by solvde.*/
{
  int loff,l,j,ic,i;
  double vx;
  loff=jc1-jm1;
  ic=ic1;
  for (j=jz1;j<=jz2;j++) { //Loop over columns to be zeroed.
    for (l=jm1;l<=jm2;l++) { //Loop over columns altered.
      vx=c[ic][l+loff][kc];
      for (i=iz1;i<=iz2;i++) s[i][l] -= s[i][j]*vx; //Loop over rows.
    }
    vx=c[ic][jcf][kc];
    for (i=iz1;i<=iz2;i++) s[i][jmf] -= s[i][j]*vx; //Plus final element.
    ic += 1;
  }
}

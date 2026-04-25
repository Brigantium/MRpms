#define R_NO_REMAP
#include "header.h"

double paso(double x, double y, double h1, double h2, double *X, double *Y, int n) {
  
  double aux = 0.;
  double res, sumaux;
  res = 0.;
  sumaux = 0.;
  for(int i = 0;i<n;i++){
    aux = exp(-pow(x-X[i],2)/(h1*h1*2.)-pow(y - Y[i],2)/(2.*h2*h2));
    res = res + aux*Y[i];
    sumaux = sumaux + aux;
  }
  res = res/ sumaux;
  return res;
}
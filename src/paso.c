#define R_NO_REMAP
#include "header.h"

double paso(double x, double y, double h1, double h2, double *X, double *Y, int n) {
  
  double aux = 0.;
  double res, sumaux;
  res = 0.;
  sumaux = 0.;
  double h1aux = h1*h1*2.;
  double h2aux = 2.*h2*h2;
  for(int i = 0;i<n;i++){
    aux = exp(-pow(x-X[i],2)/(h1aux)-pow(y - Y[i],2)/(h2aux));
    res = res + aux*Y[i];
    sumaux = sumaux + aux;
  }
  res = res/ sumaux;
  return res;
}



double pasoEf(double *xaux, double y, double h2aux, double *Y, int n) {
  
  double aux = 0.;
  double res, sumaux;
  res = 0.;
  sumaux = 0.;
  for(int i = 0;i<n;i++){
    aux = xaux[i]*exp(-pow(y - Y[i],2)/(h2aux));
    res = res + aux*Y[i];
    sumaux = sumaux + aux;
  }
  res = res/ sumaux;
  return res;
}

double *pasoauxx(double x, double *X, double h1aux, int n){

  double *res = (double *)malloc(sizeof(double) *n);

  for(int i = 0; i < n; i++){
    res[i] = exp(-pow(x-X[i],2)/(h1aux));
  }


  return res;

}


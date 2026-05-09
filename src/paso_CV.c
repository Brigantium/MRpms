#define R_NO_REMAP
#include "header.h"

double paso_CV(double x, double y, double h1, double h2, double *X, double *Y, int n, int I) {
  
  double aux = 0.;
  double res, sumaux;
  res = 0.;
  sumaux = 0.;
  // aux[I] = 0.;
  for(int i = 0; i<n;i++){
    aux = exp(-pow(x-X[i],2)/(h1*h1*2.)-pow(y - Y[i],2)/(2.*h2*h2));
    res = res + aux*Y[i];
    sumaux = sumaux + aux;
    if (i+1 == I){
      i++;
    }
  }
  res = res/ sumaux;  
  return res;
}


double paso_CVEf(double *xaux, double y, double h2aux, double *Y, int n, int I) {
  
  double aux = 0.;
  double res, sumaux;
  res = 0.;
  sumaux = 0.;
  // aux[I] = 0.;
  for(int i = 0; i<n;i++){
    aux = xaux[i]*exp(-pow(y - Y[i],2)/(h2aux));
    res = res + aux*Y[i];
    sumaux = sumaux + aux;
    if (i+1 == I){
      i++;
    }
  }
  res = res/ sumaux;  
  return res;
}



double *paso_CVauxx(double x, double h1aux, double *X, int n, int I) {
  
  double *res = (double*)malloc(sizeof(double) *n);
  res[I] = 0.;

  // aux[I] = 0.;
  for(int i = 0; i<n;i++){
    res[i] = exp(-pow(x - X[i],2)/(h1aux));

    if (i+1 == I){
      i++;
    }
  }

  return res;
}
#define R_NO_REMAP
#include "header.h"
#include "paso.h"


double *PMS1c(double *X, double *Y, double x, double *ymalla, double h1, double h2, double eps, int k, int n){

  double *y = (double*)malloc(sizeof(double)*k);
  short int yend = 0;
  short int yend_index = 0;
  double cambio_aux;

  while(1){
    for(int i = yend_index; i < k; i++){
      y[i] = paso(x,ymalla[i],h1,h2,X,Y,n);

      if (fabs(y[i] - ymalla[i]) < eps){
        cambio_aux = y[yend];
        y[yend] = y[i];
        y[i] = cambio_aux;
        yend++;
      }

    }

    if (yend > (double)k-1.5){
      break;
    }
    yend_index = yend;
    for(int i = yend_index; i<k;i++){
      ymalla[i] = y[i];
    }


  }
  return y;
}
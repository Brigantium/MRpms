#define R_NO_REMAP
#include "header.h"
#include "paso_CV.h"


double *PMS1_CV(double *X, double *Y, double x, double *ymalla, double h1, double h2, double eps, int k, int n, int I){

  double *y = (double*)malloc(sizeof(double)*k);
  short int yend = 0;
  short int yend_index = 0;
  double cambio_aux;

  while(1){
    for(int i = yend_index; i < k; i++){
      y[i] = paso_CV(x,ymalla[i],h1,h2,X,Y,n,I);

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
  return(y);

}
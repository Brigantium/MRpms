#define R_NO_REMAP
#include "header.h"
#include "paso.h"

SEXP PMS1(SEXP Xarg, SEXP Yarg, SEXP xarg, SEXP ymallaarg, SEXP h1arg, SEXP h2arg, SEXP epsarg, SEXP karg, SEXP narg){
  double *X = REAL(Xarg);
  double *Y = REAL(Yarg);
  double x = REAL(xarg)[0];
  double *ymalla = REAL(ymallaarg);
  double h1 = REAL(h1arg)[0];
  double h2 = REAL(h2arg)[0];
  double eps = REAL(epsarg)[0];
  int k = Rf_asInteger(karg);
  int n = Rf_asInteger(narg);
  
  double *y = (double*)malloc(sizeof(double)*k);
  int yend = 0;
  
  while(1){
    y[0] = paso(x,ymalla[0],h1,h2,X,Y,n);
    
    yend = (fabs(y[0] - ymalla[0]) < eps);
    
    for(int i = 1;i <k; i++){
      y[i] = paso(x,ymalla[i],h1,h2,X,Y,n);
      
      yend += (fabs(y[i] - ymalla[i]) < eps);
    }
    
    if(yend > (double)k-0.5){
      break;
    }
    
    for(int i = 0;i < k;i++){
      ymalla[i] = y[i];
    }
    yend = 0;
    
  }
  
  SEXP  yres = PROTECT(Rf_allocVector(REALSXP, k));
  for(int i = 0;i < k;i++){
    REAL(yres)[i] = y[i];
    // SET_VECTOR_ELT(yres, i, ScalarReal(y[i]));
  }
  UNPROTECT(1);
  
  free(y);

  return yres;
}
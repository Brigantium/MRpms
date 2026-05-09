#ifndef _PASOCV_
#define _PASOCV_

double paso_CV(double x, double y, double h1, double h2, double *X, double *Y, int n, int I);
double paso_CVEf(double *xaux, double y, double h2aux, double *Y, int n, int I);
double *paso_CVauxx(double x, double h1aux, double *X, int n, int I);
#endif
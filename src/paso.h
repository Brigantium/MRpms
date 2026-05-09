#ifndef _PASO_
#define _PASO_

double paso(double x, double y, double h1, double h2, double *X, double *Y, int n);
double pasoEf(double *xaux, double y, double h2aux, double *Y, int n);
double *pasoauxx(double x, double *X, double h1aux, int n);

#endif
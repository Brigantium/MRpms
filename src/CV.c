#define R_NO_REMAP
#include "header.h"
#include "PMS1_CV.h"


SEXP CV(SEXP Xarg, SEXP Yarg, SEXP mallarg, SEXP dimarg, SEXP h1arg, SEXP h2arg, SEXP parg, SEXP epsarg, SEXP narg, SEXP karg){
  double *X = REAL(Xarg);
  double *Y = REAL(Yarg);
  short int dim = Rf_asInteger(dimarg);
  double *malla = REAL(mallarg); // tendremos que acceder como un vector obtenido tras concatenar las columnas de la matriz original
  double h1 = REAL(h1arg)[0];
  double h2 = REAL(h2arg)[0];
  short int p = Rf_asInteger(parg);
  double eps = REAL(epsarg)[0];
  int k = Rf_asInteger(karg);
  int n = Rf_asInteger(narg);

  // resultado a devolver
  double rescv = 0.;

  // matriz de modas, en las filas estarán las distintas modas para cada
  //  x de la muestra
  
  double *ymalla_aux = (double*)malloc(sizeof(double)*k);
  
  short int nmodas = 0;
  
  double aux2 = 0., aux = 0.;
  
  bool moda_distinta = 0;

  double h1aux = h1*h1*2.;
  double h2aux = h2*h2*2.;
  
  
  for(int i = 0; i < n;i++){
    double *modas_aux = (double*)malloc(sizeof(double)*k);
    // inicializamos el contador de modas para el i-ésimo x:
    nmodas = 0;

    // copiamos los valores de la malla en Y que necesitamos en un puntero auxiliar
    for(int j =0; j<k;j++){
      ymalla_aux[j] = malla[i*k+j +dim*k*n];
      // Rprintf("ymalla[%i]: %f\n",j,ymalla_aux[j]);
    }

    // calculamos las modas 
    // Rprintf("x: %f --- X[i]: %f\n",malla[i*k], X[i]);
    modas_aux = PMS1_CV(X,Y,malla[i*k],ymalla_aux,h1aux,h2aux,eps, k, n, i);
    // Rprintf("modas_aux: ");
    // for (int j = 0; j <k; j++){
    //   Rprintf("%f, ",modas_aux[j]);

    // }
    // Rprintf("\n");
    // Rprintf("modas_aux: %f, %f, %f, %f, %f ",modas_aux[0],modas_aux[1],modas_aux[2],modas_aux[3],modas_aux[4]);
    
    // guardamos las modas únicas 
    double *unique_modas_aux = (double*)malloc(sizeof(double)*k);

    // comenzamos a limpiar los valores para quedarnos con las modas procesadas
    unique_modas_aux[nmodas] = (double)round(modas_aux[0]*(double)pow(10.,(p-2)))/(double)pow(10.,(p-2));
    nmodas++; //contamos la primera moda
    
    
    for(int j = 1; j<k;j++){
      modas_aux[j] = (double)round(modas_aux[j]*(double)pow(10.,(p-2)))/(double)pow(10.,(p-2)); // procesamos la moda
      
      moda_distinta = 1;
      for(int it = 0;it<nmodas;it++){
        if (fabs(modas_aux[j] - unique_modas_aux[it])<eps){
          moda_distinta = 0;
          break;
        } 
        
      }
      if(moda_distinta){
        unique_modas_aux[nmodas] = modas_aux[j];
        nmodas++;
        
      }
      
    }
    
    // Rprintf("nmodas %i\n",nmodas);
    // procedemos con el cálcumo del error:
    aux2 = INFINITY;
    for (int j = 0; j<nmodas; j++){
      aux = pow(Y[i] - unique_modas_aux[j],2);

      // hacemos la búsqueda del mínimo de las distancias
      if (aux < aux2){
        aux2 = aux;
      }
    }

    // actualizamos el resultado final
    rescv = rescv + (double)(nmodas*nmodas)*aux2;
    free(modas_aux);
    free(unique_modas_aux);
  }
  // liberamos la memoria reservada
  
  free(ymalla_aux);

  rescv = rescv/(double)n;

  SEXP res = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(res)[0] = rescv;
  UNPROTECT(1);
  return(res);
}
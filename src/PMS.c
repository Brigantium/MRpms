#define R_NO_REMAP
#include "header.h"
#include "PMS1c.h"

SEXP PMS(SEXP Xarg, SEXP Yarg, SEXP mallarg, SEXP dimarg, SEXP h1arg, SEXP h2arg, SEXP parg, SEXP epsarg, SEXP narg, SEXP karg, SEXP larg){

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
  int l = Rf_asInteger(larg);

  SEXP modasres = PROTECT(Rf_allocVector(VECSXP, l));


  double *ymalla_aux = (double*)malloc(sizeof(double)*k);
  
  short int nmodas = 0;
  
  bool moda_distinta = 0;


  for(int i = 0; i <l;i++){

    double *modas_aux = (double*)malloc(sizeof(double)*k);
    // inicializamos el contador de modas para el i-ésimo x:
    nmodas = 0;

    // copiamos los valores de la malla en Y que necesitamos en un puntero auxiliar
    for(int j =0; j<k;j++){
      ymalla_aux[j] = malla[i*k+j +dim*k*l];
      
    }

    // calculamos las modas 
    modas_aux = PMS1c(X,Y,malla[i*k],ymalla_aux,h1,h2,eps, k, n);
    

    // guardamos las modas únicas en un nuevo vector
    double *unique_modas_aux = (double*)malloc(sizeof(double)*k);

    // comenzamos a limpiar los valores para quedarnos con las modas procesadas
    unique_modas_aux[nmodas] = (double)round(modas_aux[0]*(double)pow(10.,(p-2)))/(double)pow(10.,(p-2));
    nmodas++; //contamos la primera moda
    
  
    // realizamos le proceso de búsqueda de las modas únicas
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

    
    SEXP modasi_aux = Rf_allocVector(REALSXP, nmodas);
    for(int it = 0; it <nmodas; it++){
      REAL(modasi_aux)[it] = unique_modas_aux[it];
    }
    
    SET_VECTOR_ELT(modasres, i, modasi_aux);

    free(modas_aux);
    free(unique_modas_aux);

  }


  free(ymalla_aux);
  UNPROTECT(1);
  return modasres;
}
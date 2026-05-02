#define R_NO_REMAP
#include "header.h"
#include "mtwister.h"
#include "PMS1c.h"

SEXP bootstrap_deltas(SEXP Xarg, SEXP Yarg, SEXP modasarg, SEXP nmodasarg, SEXP mallaarg, SEXP narg, SEXP karg, SEXP larg, SEXP h1arg, SEXP h2arg, SEXP parg, SEXP epsarg, SEXP dimarg, SEXP Barg, SEXP seed){

    double *X = REAL(Xarg);
    double *Y = REAL(Yarg);
    double *modas = REAL(modasarg);
    int *nmodas = INTEGER(nmodasarg);
    short int dim = Rf_asInteger(dimarg);
    double *malla = REAL(mallaarg); // tendremos que acceder como un vector obtenido tras concatenar las columnas de la matriz original
    double h1 = REAL(h1arg)[0];
    double h2 = REAL(h2arg)[0];
    short int p = Rf_asInteger(parg);
    double eps = REAL(epsarg)[0];
    int k = Rf_asInteger(karg);
    int n = Rf_asInteger(narg);
    int l = Rf_asInteger(larg);
    int B = Rf_asInteger(Barg);

    MTRand r = seedRand(Rf_asInteger(seed));
//   genRand(&r);

    double *Xb = (double*)malloc(sizeof(double)*n);
    double *Yb = (double*)malloc(sizeof(double)*n);
    double *ymalla_aux = (double*)malloc(sizeof(double)*k);

    SEXP deltasbx = PROTECT(Rf_allocVector(REALSXP, l*B));

    for(int b = 0; b<B;b++){
        for(int i = 0;i<n;i++){
            double ind_aux = genRand(&r);
            int ind = floor(n*ind_aux);
            Yb[i] = Y[ind];
            Xb[i] = X[ind];
        }


  
        short int nmodasb = 0;
        
        bool moda_distinta = 0;

        // inicialización de un indicador de posición, para guardar el índice de la primera moda en x.malla[i] de los datos originales 
        int pun = 0;

        for(int i = 0; i <l;i++){


            double *modasb_aux = (double*)malloc(sizeof(double)*k);
            // inicializamos el contador de modas para el i-ésimo x:
            nmodasb = 0;

            // copiamos los valores de la malla en Y que necesitamos en un puntero auxiliar
            for(int j =0; j<k;j++){
            ymalla_aux[j] = malla[i*k+j +dim*k*l];
            
            }

            // calculamos las modas en x.malla[i]
            modasb_aux = PMS1c(Xb,Yb,malla[i*k],ymalla_aux,h1,h2,eps, k, n);
            // Rprintf("modasb_aux: ");
            // for (int auxit = 0;auxit<k;auxit++){
            //     Rprintf(" %f,",modasb_aux[auxit]);

            // }
            // Rprintf("\n");

            // guardamos las modas únicas en un nuevo vector
            double *unique_modasb_aux = (double*)malloc(sizeof(double)*k);

            // comenzamos a limpiar los valores para quedarnos con las modas procesadas
            unique_modasb_aux[nmodasb] = (double)round(modasb_aux[0]*(double)pow(10.,(p-2)))/(double)pow(10.,(p-2));
            nmodasb++; //contamos la primera moda
            
            double nuevo_minimo = INFINITY;
            double nuevo_candidato;
            for(int it = 0;it <nmodas[i];it++){
                nuevo_candidato = fabs(modasb_aux[0] - modas[it + pun]);
                // Rprintf("%f\n",nuevo_candidato);
                if (nuevo_candidato < nuevo_minimo){
                    nuevo_minimo = nuevo_candidato;
                }

            }
        
            // realizamos le proceso de búsqueda de las modas únicas, a la vez que computamos la distancia Hausdorff con las modas originales en x.malla[i]
            double nuevo_maximo = nuevo_minimo;
            for(int j = 1; j<k;j++){
                modasb_aux[j] = (double)round(modasb_aux[j]*(double)pow(10.,(p-2)))/(double)pow(10.,(p-2)); // procesamos la moda
                
                moda_distinta = 1;
                for(int it = 0;it<nmodasb;it++){
                    if (fabs(modasb_aux[j] - unique_modasb_aux[it])<eps){
                    moda_distinta = 0;
                    break;
                    } 
                    
                }
                if(moda_distinta){
                    unique_modasb_aux[nmodasb] = modasb_aux[j];
                    nmodasb++;

                    // ahora buscamos el mínimo de las diferencias de esta nueva moda con todas las modas de los datos originales
                    double nuevo_minimo = INFINITY;
                    double nuevo_candidato;
                    for(int it = 0;it <nmodas[i];it++){
                        nuevo_candidato = fabs(modasb_aux[j] - modas[it + pun]);
                        // Rprintf("%f\n",nuevo_candidato);
                        if (nuevo_candidato < nuevo_minimo){
                            nuevo_minimo = nuevo_candidato;
                        }

                    }


                    // y computamos la distancia Hausdorff buscando el máximo de las mínimas distancias entre las modas 
                    if(nuevo_minimo >nuevo_maximo){
                        nuevo_maximo = nuevo_minimo;
                    }
                }
            
            }
            pun += nmodas[i];
            free(unique_modasb_aux);
            free(modasb_aux);

            REAL(deltasbx)[l*b + i] = nuevo_maximo;
        }

    }
    free(ymalla_aux);
    free(Xb);
    free(Yb);

    UNPROTECT(1);
    // deltas va a ser un único gran vector de longitud l*B
    return deltasbx;
}

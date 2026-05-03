#' Construct a mesh compatible with MRpms functions.
#' 
#' Given a matrix with data points, or two matrix with the covariable points and response values,
#'  construct a mesh of points in both ways: compute ´l´ diferent equidistant `x` points and asign `k` diferent `y`
#'  values to each ones between the nearest data points; or take a given `x` matrix of points and preforms the same construction.
#' 
#' @param muestra Matrix with data points.
#' 
#' @param X Matrix with covariable data points.
#' 
#' @param Y Vector with response variable data values.
#' 
#' @param k Number of `y` values for each `x` point in the new mesh.
#' 
#' @param len Number of diferent `x` points in the new mesh.
#' 
#' @param dim Dimension of covariable.
#' 
#' @param x.malla Matrix of `x` points given by user. 
#' 
#' @return A matrix with the new mesh. Each `x` point is repeated `k` times, each one with 
#' their associated `y` value. Besides, the matrix have three attributes: `l`, `k` and `x.malla`
#'  which contains a matrix with the diferents `x` points.
#' @export


mallador <- function(X = muestra[,1:dim], Y = muestra[,dim+1], 
                    k = 10, len = 100, muestra = cbind(X,Y), dim = 1, x.malla){

  if(dim > 1){ # caso multidimensional
    
    if(missing(x.malla)){ # en caso de no recibir unos puntos predefinidos, construimos la malla de la covariable mediante puntos equidistantes
      minx <- apply(X, 2, min)
      maxx <- apply(X, 2, max)
      x.temp <- apply(cbind(minx, maxx), 1, function(x) seq(x[1], x[2], len = len))
      x.malla <- sapply(1:dim, function(i) rep(rep(x.temp[,i], each = len^(dim-i)), times = len^(i-1)))
    } else {
      len <- nrow(x.malla)
    }
    
    y.temp <- c(apply(x.malla, 1 ,function(x) Y[order(MRpms_norma(x, X))[1:k]])) # buscamos los valores de Y de los puntos de la muestra más cercanos a `x`
    malla <- cbind(matrix(rep(x.malla, each = k), ncol = dim), y.temp) # construimos la malla

  }else{ # caso unidimensional
    
    if(missing(x.malla)){ # en caso de que no recibir unso puntos predefinidos, construimos la malla de la covariable mediante putnos equidistantes
      minx <- min(X)
      maxx <- max(X)
      x.malla <- seq(minx, maxx, len = len)
    } else { # en otro caso, 

      len <- length(x.malla)
    }
    
    X.aux <- matrix(X,ncol=1) # reescribimos X para que pueda ser usada en la función auxiliar norma
    
    y.temp <- c(sapply(x.malla, function(x) Y[order(MRpms_norma(x, X.aux))[1:k]]))
    x.temp <- rep(x.malla, each = k)
    malla <- cbind(x.temp,y.temp)
  }
  

  malla <- structure(malla, 
                      len = len,
                      k = k,
                      x.malla = x.malla,
                      class = "MRpms_malla")
  return(malla = malla)
}
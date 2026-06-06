#' @title Constructs a mesh compatible with MRpms functions
#'
#' @description Given a matrix containing the sample, constructs a mesh.
#' This can be done in two ways: either computing `l` different equidistant `x`
#' points and assigning `k` different `y` values to each one,
#' chosen in a k nearest neighbors fashion or, if given a matrix of `x`
#' points (`x.malla`), performing the same construction over them.
#'
#' @param data Matrix or data frame containing the sample.
#'
#' @param k Number of `y` values for each `x` point in the new mesh.
#'
#' @param len Number of different `x` points in the new mesh.
#'
#' @param x.malla Matrix of `x` points given by user.
#'
#' @return A matrix containing the points of the new mesh.
#' Each `x` point is repeated `k` times, each one with an associated `y` value.
#' This matrix has three attributes: `len`, `k` and `x.malla`, the last of which
#' is a matrix containing the different `x` points in the mesh.
#'
#' @examples
#' malla <- mallador(twosines)
#' plot(twosines)
#' points(malla, pch = 19, cex = 0.8)
#'
#' @export


mallador <- function(data, k = 10, len = 200, x.malla){

  dim = ncol(data) - 1
  X = data[, 1:dim]
  Y = data[, dim+1]

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



MRpms_norma <- function(v,V) sqrt(apply((v-V)^2,1,sum))
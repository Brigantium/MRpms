#' @title Compute the volume of the uniform prediction set.
#' 
#' @description Compute the volume of the uniform prediction set for a fitted model. First compute the modes in a given mesh 
#' using cross-validation, compute the quantile for a given confidence level and later use Monte-Carlo integration. Currently,
#' only works in the univariate case.
#' 
#' @param X Matrix with the coordinates of the covariable vector.
#' 
#' @param Y Vector with the response variable values.
#' 
#' @param malla A `malla` object from the `mallador` function.
#' 
#' @param dim Dimension of covariable vector.
#' 
#' @param h1 Smooth parameter for covariable.
#' 
#' @param h2 Smooth parameter for response variable.
#' 
#' @param p Number of decimal precision.
#' 
#' @param eps Tolerance for convergence in the partial Mean-Shift.
#' 
#' @param n Sample size.
#' 
#' @param conf.level Confidence level for the uniform prediction set.
#' 
#' @param k Number of initial guesses in every covariable point in the mesh.
#' 
#' @param lensopx Size of the support of covariable.
#' 
#' @return The volume of the estimated uniform prediction set.
#' 
#' @keywords internal


hP <-function(X, Y, malla, dim = 1, h1, h2,
              p = floor(-log(eps, base = 10)), eps = 10^(-p),
              n = length(Y), conf.level = 0.95, k = attr(malla, "k"),
              lensopx = max(X) - min(X)){

  # para cada punto de la muestra calculamos sus modas sin tenerlo en cuenta.
  aux <- lapply(1:n, function(i) sort(unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
                                                                x = malla[(i-1)*k+1,1:dim],
                                                                ymalla = malla[(i-1)*k+(1:k),dim+1],
                                                                h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-2))))
  
  # tomamos la diferencia entre las modas que abarca el nivel de confianza pedido de los valores de la respuesta en todos los puntos
  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), conf.level)

  # realizamos la integración monte-carlo.
  ni <- sapply(aux, length)
  intk <- mean(ni)*lensopx

  return(epsh*intk)
}
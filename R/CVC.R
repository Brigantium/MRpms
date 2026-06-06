#' @title Compute the cross-validation like function suggested by Zhou and Huang
#' 
#' @description This function mask a C function which preforms the cross validation 
#' suggested by Zhou and Huang to search the optimal smooth parameter.
#' 
#' @param data Matrix or data frame containing the sample.
#' 
#' @param X Matrix with the covariable points.
#' 
#' @param Y Vector with the response values.
#' 
#' @param malla An `MRpms_malla` object returned by the `mallador` function with the initial guesses to comppute the modes.
#' 
#' @param dim Dimension of covariable.
#' 
#' @param h1 Smooth parameter for covariable.
#' 
#' @param h2 Smooth parameter for response.
#' 
#' @param p Decimal precision.
#' 
#' @param eps Convergence tolerance.
#' 
#' @param n Sample size.
#' 
#' @param k Number of guesses for each covariable point in mesh.
#' 

#' @returns The value of the cross-validation like function suggested by Zhou and Huang,
#' 
#' \deqn{\displaystyle \frac{1}{n}\sum_{i=1}^n d^2\left(\widehat{M}_{n,\mathbf{h},-i}\left(X_i\right),Y_i\right)\widehat{N}_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right),}
#' 
#' where \eqn{\widehat M_{n,\mathbf{h},-i}(X_i)} is the estimated modes in \eqn{X_i} without consider \eqn{X_i} and \eqn{\widehat{N}_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right)}
#' is the number of modes estimated in \eqn{X_i}.
#' 
#' 
#' @references 
#' Zhou, H. and Huang, X. (2019).
#' *Bandwidth selection for nonparametric modal regression*.
#' Communication in Statistics - Simulation and Computation,
#' **48**(4), 968--984.
#' 
#' @keywords internal  
## usethis namespace: start
#' @useDynLib MRpms, .registration = TRUE
## usethis namespace: end

CVC <- function(data = cbind(X,Y),X = data[,1:dim],Y = data[,dim+1],
                malla = mallador(data), dim, h1, h2, p = floor(-log(eps, base = 10)),
                eps = 10^(-p), n = length(Y), k = attr(malla,"k")){

  .Call(CV, as.numeric(X), as.numeric(Y), as.numeric(malla),
              as.integer(dim),as.numeric(h1),as.numeric(h2),
              as.integer(p),as.numeric(eps),as.integer(n), as.integer(k))
}

#' Pointwise and uniform confidence intervals.
#' 
#' Compute pointwise and uniform confidence intervals to each mode. Besides, plot both confidence intervales and data points.
#' 
#' @param muestra Matrix which contains the data points.
#' 
#' @param X Matrix of covariable data points in case `muestra` not present.
#' 
#' @param Y Vector with response values in case `muestra` not present.
#' 
#' @param dim DImension of covariable.
#' 
#' @param modas List of vector which contains the modes on each covariable point of `malla`. In case `modas` missing, the function calls `PMSc` to compute them.
#' 
#' @param malla Matrix of points used to compute the modes using the Partial Mean-Shift algorithm. If not provided, function calls  `mallador` to compute 
#' 
#' @param h1 Smoothing parameter for covariables.
#' 
#' @param h2 Smoothing parameter for response variable.
#' 
#' @param type Numerical argument. "0" to compute pointwise intervals, "1" to uniform intervals or "2" to compute both. By default, preforms both.
#' 
#' @param eps Tolerance of convergence and minimum diference between modes.
#' 
#' @param p Number of significance digits considered to diference modes.
#' 
#' @param conf.level Confidence level of intervals.
#' 
#' @param B Number of _bootstrap_ simulations to compute the confidence intervals.
#' 
#' @param seed Seed to generate the _bootstrap_ simulations.
#' 
#' @return A list. In all three cases of `type` values firts entry are a point list with all modes. If `type` is 0, returns a vector with every pointwise ratios too; if `type` is 1 returns the uniform ratio.
#' 
#'
#' @export


ConfPMS <- function(muestra = cbind(X,Y), modas, dim = ncol(muestra) -1, 
                    malla = mallador(X, Y, dim = dim), h1 = 0.3, h2 = 0.5, 
                    eps = 1e-8, p = floor(-log(eps, base = 10)),
                    X = muestra[,1:dim], Y = muestra[,dim+1], 
                    conf.level=0.95, B=500, type = 2,
                    seed = 2026){

  
  x.malla <- attr(malla,"x.malla")
  k <- attr(malla, "k")
  l <- attr(malla, "len")
  n <- length(Y)
  
  if(missing(modas)){

    modas <- PMSc(X,Y,malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                  n = n, k = k, l = l)
  }
  plot(X, Y)
  ni <- sapply(modas, length)
  xg <- rep(x.malla, times = ni)
  yg <- unlist(modas)
  salida <- list(puntos=cbind(xg,yg))
  
  X <- matrix(X, ncol=dim)
  set.seed(seed)
  Deltasbx <- matrix(Bdeltas(X,Y,modas,ni,malla,n,k,l,h1,h2,p,eps,dim,B,seed),nrow = l, byrow = TRUE)
  # Deltasbx <- replicate(B, .Deltas(X = X, Y = Y, modas = modas,
  #                                  malla = malla, n = n, k = k, l = l, 
  #                                  h1 = h1, h2 = h2, p = p, eps = eps, 
  #                                  dim = dim))
  
  if(type==1 | type==2){
    Deltasb <- apply(Deltasbx, 2, max)
    delta <- stats::quantile(Deltasb, conf.level)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-delta, yg[i]+delta), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - delta, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + delta, col = "grey", pch = 19, cex = 0.4)
    salida <- append(salida, list(delta = delta))
  }
  
  if(type==0 | type==2){
    deltasx <- apply(Deltasbx, 1, function(x) stats::quantile(x, conf.level))
    radio <- rep(deltasx, ni)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-radio[i], yg[i]+radio[i]), col ="blue", lwd = 1.5))
    salida <- append(salida, list(deltax = radio))
  }
  
  graphics::points(xg, yg, col = "red", pch = 19)
  
  return(salida)
}
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
#' @param modas List of vector which contains the modes on each covariable point of `malla`. In case `modas` missing, the function calls `PMSc` to compute them. It must be an special class object MRpms_modas.
#' 
#' @param malla Matrix of points used to compute the modes using the Partial Mean-Shift algorithm. If not provided, function calls  `mallador` to compute. It must be an special class object MRpms_malla.
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
#' @param k Number of Y values per x value in the `malla` object created, if `malla` is not provided.
#' 
#' @param l Number of X diferent points computed in the `malla` object.
#' 
#' @return A list. In all three cases of `type` values, firts entry is a list with all modes estimated. If `type` is 0,
#'  returns a vector with every pointwise ratios too; if `type` is 1 returns the uniform ratio.
#' 
#'
#' @export


ConfPMS <- function(muestra = cbind(X,Y), modas,  
                    malla, h1 = 0.3, h2 = 0.5, 
                    eps = 1e-8, p = floor(-log(eps, base = 10)),
                    k , l,
                    conf.level=0.95, B=500, type = 2,
                    seed = 2026){
                      
  
  # comprobación de los argumentos suministrados:
  if(!methods::is(h1,"numeric") | h2==0){
    stop("Bandwidth for covariable estimator `h1` is not a positive (>0) numeric value.")
  }

  
  if(!methods::is(h2,"numeric")| h2==0){
    stop("Bandwidth for response estimator `h2` is not a positive (>0) numeric value.")
    
  }
  
  if(!methods::is(p,"numeric")){
    stop("Precision must be an unidimensional numeric value.")
  }
  p = floor(abs(p))
  
  if(!methods::is(k,"numeric")){
    stop("The number of initial Y values per x.malla point must be an unidimensional integer value.")
  }
  k = floor(abs(k))
  
  
  if (!methods::is(conf.level, "numeric") | length(as.vector(conf.level)) != 1L | conf.level > 1 | conf.level <0){
    stop("Confidence level must be a number between 0 and 1.")
  }
  
  

  if(!missing(malla)){
    if(!methods::is(malla, "MRpms_malla")) stop("`malla` is not a `malla` object. Use `mallador` function to compute one.")
  }

  # comprobar si recibimos una muestra válida

  # comprobamos si de hecho es una matriz o array
  if((!methods::is(muestra,"array") & !methods::is(muestra,"matrix")) | typeof(muestra) != "double" ){
    stop("`muestra` argument is not a numeric matrix nor an array.")
    
  }
  
  
  # calculamos la dimensión de la covariable
  dim = ncol(muestra) - 1

  # comprobamos que el número de dimensiones sea correcto
  if (!methods::is(dim,"numeric")){
    stop("Number of dimensions not correct. Pleas, check `muestra` has more than two columns.")
    
  }

  # Separamos la variable explicativa de la variable respuesta
  X = muestra[,1:dim]
  Y = muestra[,dim+1]

  # si no se da una malla, se especifica una  
  if(missing(malla)){
    
    # de ser proveído un objeto modas, entonces construimos la malla con los puntos con los que esta fue calculada
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")
      if (missing(k)){
        stop("Provide the desired number of Y values per x point on the `malla`, `k`.")
        
      }
      malla = mallador(X,Y, dim = dim, x.malla = attr(modas,"x.malla"), k = k)
    } else{  # si no fue provisto un objeto `modas`, usamos los argumentos suministrados
      if(missing(k) | missing(l)){
        stop("Not enough arguments to compute a `malla` object. Please, check `k` and `l` argument were provided.")
        
      }
      malla = mallador(X, Y, dim = dim, k = k, len = l)
    }
  } else {
    # en caso de que tengamos ambos, comprobemos que están definidos sobre los mismos 
    # puntos.
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")
      if(attr(modas,"x.malla") != attr(malla,"x.malla")){
        stop("`x.malla` in both `malla` object and `modas` object are not equals.
Please, provide a `malla` object build over the same `x.malla` as `modas`.")
        
      }
    }
  }
  
  # obtenemos la información de la malla necesaria
  x.malla <- attr(malla,"x.malla") # puntos sobre los que está construida
  k <- attr(malla, "k") # número de Y's por punto de la malla
  l <- attr(malla, "len") # número de puntos x diferentes en la malla
  n <- length(Y) # tamaño muestral
  

  if(missing(modas)){

    modas <- PMSc(X,Y,malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                  n = n, k = k, l = l)
  }

  modas <- structure(modas, 
                    x.malla = x.malla,
                  class = "MRpms_modas")
  # graficamos las modas calculadas
  plot(X, Y)
  ni <- lapply(modas, length)
  xg <- rep(x.malla, times = ni)
  yg <- unlist(modas)
  # salida <- list(puntos=cbind(xg,yg))
  salida <- list(modas = modas, x.malla = x.malla)

  X <- matrix(X, ncol=dim)
  set.seed(seed)
  # Deltasbx <- matrix(Bdeltas(X,Y,modas,ni,malla,n,k,l,h1,h2,p,eps,dim,B,seed),nrow = l, byrow = TRUE)
  Deltasbx <- replicate(B, .Deltas(X = X, Y = Y, modas = modas, 
                                   malla = malla, n = n, k = k, l = l, 
                                   h1 = h1, h2 = h2, p = p, eps = eps, 
                                   dim = dim))
  
  if(type==1 | type==2){
    Deltasb <- apply(Deltasbx, 2, max)
    delta <- stats::quantile(Deltasb, conf.level)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-delta, yg[i]+delta), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - delta, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + delta, col = "grey", pch = 19, cex = 0.4)
    salida <- append(salida, list(delta = delta))
  } else if(type==0 | type==2){
    deltasx <- apply(Deltasbx, 1, function(x) stats::quantile(x, conf.level))
    radio <- rep(deltasx, ni)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-radio[i], yg[i]+radio[i]), col ="blue", lwd = 1.5))
    salida <- append(salida, list(deltax = radio))
  } else {
    stop("Please, select `type` between 0, 1 or 2.")
  }

  # graphics::plot(modas, pch = 19, col = "red")
  graphics::points(xg, yg, col = "red", pch = 19)
  
  return(salida)
}
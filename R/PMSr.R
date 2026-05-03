#' Function to performs modal regression
#' 
#' This function performs modal regression over a given sample of unidimensional, in future multidimensional, covariable, using a 
#' given set of points where estimate.
#' 
#' @param muestra Matrix which contains the data points. 
#' 
#' @param x.malla Set of covariable points where compute the modes. 
#' 
#' @param h1 Smoothing parameter for covariables.
#' 
#' @param h2 Smoothing parameter for response variable.
#' 
#' @param eps Tolerance of convergence and minimum diference between modes.
#' 
#' @param p Number of significance digits considered to diference modes.
#' 
#' @param k Number of Y values per x value in the `malla` object created, if `malla` is not provided.
#' 
#' @param l Number of X diferent points computed in the `malla` object.
#' 
#' @param malla Matrix of points used to compute the modes using the Partial Mean-Shift algorithm. If not provided, function calls  `mallador` to compute.
#'
#' @return The function returns a list of `l` entries, each of one contains the modes estimated in each `x.malla`. In addition, `x.malla` is returned as attribute.
#' 
#' @export

PMS <- function(muestra, x.malla,
                h1 = 0.5, h2 = 1, 
                p = floor(-log(eps,base = 10)), eps = 10^(-p), k ,l,  malla){
  

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


  # comprobamos si de hecho es una matriz o array
  if((!methods::is(muestra,"array") & !methods::is(muestra,"matrix")) | typeof(muestra) != "double" ){
    stop("`muestra` argument is not a numeric matrix nor an array.")
    
  }

  # Comprobamos si x.malla es un vector/matriz de puntos numérico
  if (!missing(x.malla)){
    if((!methods::is(x.malla,"array") & !methods::is(x.malla,"matrix") & !methods::is(x.malla,"numeric")) | typeof(x.malla) != "double" ){
      stop("`x.malla` argument is not a numeric matrix nor an array.")
    }
  }
  
  dim = ncol(muestra) - 1
  
  if(!methods::is(dim,"numeric")){
    stop("Incorrect dimension. Please, make sure `muestra` is a matrix with at least two columns, one with te X values and the other 
    with the Y ones.")
    
  }
  
  
  # en caso de que no sea dado un objeto malla, se calcula
  if (missing(malla)){
    if(missing(k)){
      stop("Please, provide the desired number of Y values per x point in the `malla` object, `k`.")
      
    } 
    if (missing(l) & missing(x.malla)){
      stop("Please, provided the desired number of diferent x point in the `malla` object where the modes will be estimated, `l`.
      As an alternative to `l`, you can provide an array with the desired x points with the `x.malla` argument.")
      
    }
    if (!missing(l)) malla = mallador(muestra[,1:dim], muestra[,dim+1],dim = dim, k = k, len = l)
    if (!missing(x.malla)) malla = mallador(muestra[,1:dim], muestra[,dim+1],dim = dim, k = k, x.malla = x.malla)
      
  } 
  
  # paso = if(dim==1) .pasoU else .pasoM
  l = attr(malla, "len")
  k = attr(malla,"k")
  n = nrow(muestra)

  # realizamos la estimación de las modas
  modas <- PMSc(muestra[,1:dim], muestra[,dim+1], malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                           n = n, k = k, l = l)


  modas <- structure(modas, 
                    x.malla = attr(malla,"x.malla"),
                  class = "MRpms_modas")
  return(modas)
}

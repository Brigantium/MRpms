#' Performs modal regression over a given sample of unidimensional or multidimensional covariable.
#' 
#' @param muestra Matrix which contains the data points. 
#' 
#' @param malla Matrix of points used to compute the modes using the Partial Mean-Shift algorithm. If not provided, function calls  `mallador` to compute.
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
#' @param x.malla Set of covariable points where compute the modes. 
#' 
#' @export

PMS <- function(muestra, x.malla,
                h1 = 0.5, h2 = 1, 
                p = floor(-log(eps,base = 10)), eps = 10^(-p), k ,l,  malla){
  
  
  dim = ncol(muestra) - 1
  
  if(!methods::is(dim,"numeric")){
    message("Incorrect dimension. Please, make sure `muestra` is a matrix with at least two columns, one with te X values and the other 
    with the Y ones.")
    return()
  }
  

  # en caso de que no sea dado un objeto malla, se calcula
  if (missing(malla)){
    if(missing(k)){
      warning("Please, provide the desired number of Y values per x point in the `malla` object, `k`.")
      return()
    } 
    if (missing(l) & missing(x.malla)){
      warning("Please, provided the desired number of diferent x point in the `malla` object where the modes will be estimated, `l`.
       As an alternative to `l`, you can provide an array with the desired x points with the `x.malla` argument.")
      return()
    }
    if (!missing(l)) malla = mallador(muestra[,1:dim], muestra[,dim+1],dim = dim, k = k, len = l)
    if (!missing(x.malla)) malla = mallador(muestra[,1:dim], muestra[,dim+1],dim = dim, k = k, x.malla = x.malla)
  } 
  
  # paso = if(dim==1) .pasoU else .pasoM
  l = attr(malla, "len")
  k = attr(malla,"k")
  n = nrow(muestra)

  modas <- modas <- PMSc(muestra[,1:dim], muestra[,dim+1], malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                           n = n, k = k, l = l)


  # modas <- sapply(1:l, function(i) sort(unique(round(unlist(PMS1c(X = X, Y = Y, x = malla[(i-1)*k +1,1:dim], 
  #                                                                             ymalla = malla[(i-1)*k+(1:k),dim+1], h1 = h1, 
  #                                                                             h2 = h2,eps = 1e-8,k = 10, n = n)),p-1))))
  attr(modas,"x.malla") <- attr(malla,"x.malla")
  return(modas)
}

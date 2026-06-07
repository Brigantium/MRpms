#' Preforms the Partial Mean-Shift Algorithm
#'
#' Finds the local modes of `data` conditional on a given set of covariate values.
#'
#' @param data Matrix or data frame containing the sample. Currently, the last column is needed to be the response values.
#'
#' @param x.malla Set of covariate points where the modes will be estimated.
#'
#' @param h1 Smoothing parameter for covariates.
#'
#' @param h2 Smoothing parameter for response variable.
#'
#' @param eps Convergence tolerance.
#'
#' @param k Number of Y values per x point in the `malla` object created.
#'
#' @param len Number of different X points in the `malla` object.
#' 
#' @param dim.y Response dimension.
#'
#' @param malla Matrix containing the initial points used to compute the modes
#' with the Partial Mean-Shift algorithm. If not provided,
#' `mallador`function is called. It must be a `MRpms_malla` class object.
#'
#' @return A `MRpms_modas` object, which contains
#' * `$modas`: contains the estimated local modes conditional on the values
#'  of the covariate appearing in `x.malla`; 
#' * `$x.malla`: with each covariable point where modes were computed;
#' * `$dims`: Covariable and response dimensions. 
#'
#' @references
#' Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L. (2016).
#' *Nonparametric modal regression*. The Annals of Statistics, **44**(2), 489--514.
#'
#' @examples
#' modas <- PMS(twosines)
#' plot(modas,twosines,pch = 19, col = "red")
#'
#' @export

PMS <- function(data, x.malla = NULL,
                h1 = 0.3, h2 = 0.5, eps = 1e-8, k = 10 , len = 200, dim.y = 1,
                malla = NULL){


  if(!methods::is(h1,"numeric") | h1<=0){
    stop("Bandwidth for covariates, `h1`, is not a positive (>0) numeric value.")
  }

  if(!methods::is(h2,"numeric")| h2<=0){
    stop("Bandwidth for response variable, `h2`, is not a positive (>0) numeric value.")

  }

  if(!methods::is(eps,"numeric") | length(eps) != 1L){
    stop("Precision must be a one-dimensional numeric value.")
  }
  # p = floor(abs(p))
  p = floor(-log(eps,base = 10))

  if(!methods::is(k,"numeric")){
    stop("The number of initial Y values per x.malla point must be a one-dimensional integer value.")
  }
  k = floor(abs(k))

  if(methods::is(data,"data.frame")){ data <- as.matrix(data)}
  # comprobamos si de hecho es una matriz o array
  if((!methods::is(data,"array") & !methods::is(data,"matrix")) | typeof(data) != "double" ){
    stop("`data` argument is neither a numeric matrix nor an array.")

  }

  # Comprobamos si x.malla es un vector/matriz de puntos numérico
  if (!is.null(x.malla)){
    if((!methods::is(x.malla,"array") & !methods::is(x.malla,"matrix") & !methods::is(x.malla,"numeric")) | typeof(x.malla) != "double" ){
      stop("`x.malla` argument is neither a numeric matrix nor an array.")
    }
  }

  dim.x = ncol(data) - dim.y


  # en caso de que no sea dado un objeto malla, se calcula
  if (is.null(malla)){

    if (!is.null(x.malla)){ # si se proporcionan unos puntos donde calcular las modas, se calculan para el cálculo de la malla

      malla = mallador(data, k = k, x.malla = x.malla)

    }else malla = mallador(data, k = k, len = len)

  }

  # extraemos la información importante de la malla
  len = attr(malla, "len")
  k = attr(malla,"k")
  n = nrow(data)

  # realizamos la estimación de las modas
  modas <- PMSc(data[,1:dim.x], data[,dim.x+1:dim.y], malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim.x,
                           n = n, k = k, len = len)

  
  dims.aux <- c(dim.x,dim.y)
  names(dims.aux) <- c("X", "Y")
  hs <- c(h1,h2)
  names(hs) <- c("X","Y")
  # construimos el objeto MRpms_modas
  modas <- structure(
                    list( modas = modas,
                    x.malla = attr(malla,"x.malla"),
                    dims = dims.aux,
                    h <- hs
                  ),
                    class = "MRpms_modas")
  return(modas)
}

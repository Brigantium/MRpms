#' Uniform prediction sets based on modal regression
#'
#' Computes a uniform prediction set for the response value.
#'
#'
#' @param data Matrix or data frame containing the sample. Currently,
#' the last column is needed to be de response values.
#'
#' @param modas List containing the estimated set of conditional local
#' modes on each covariate point in `malla`. In case modas missing,
#' PMSc is called to compute them. It must be a `MRpms_modas` class object.
#'
#' @param x.malla Set of covariate points where the modes will be estimated.
#'  Not needed if `modas` is provided.
#'
#' @param h1 Smoothing parameter for covariates.
#'
#' @param h2 Smoothing parameter for response variable.
#' 
#' @param dim.y Response dimension.
#'
#' @param eps Convergence tolerance and threshold to discriminate between modes.
#'
#' @param conf.level Confidence level for the prediction set.
#'
#' @param k Number of Y values per `x` point in the `malla` object created,
#'  if `malla` is not provided.
#'
#' @param len Number of different `x` points in the `malla` object.
#'
#' @param malla Matrix containing the initial points used to compute the
#'  modes using the Partial Mean-Shift algorithm. If not provided,
#'  `mallador` function is called. It must be a MRpms_malla class object.
#'
#' @return An `MRpms_pred` object which contains,
#' * `modas`: an `MRpms_modas` object with the adjusted modes,
#' * `epsh`: the radius of the uniform prediction set.
#' * `conf.level`: the confidence level used to compute `epsh`.
#' 
#' Aditionally, the prediction set and the estimated modes are plotted in case 
#' `dim.x` and `dim.y` are `1`.
#'
#' @references
#' Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L. (2016).
#' *Nonparametric modal regression*. The Annals of Statistics, **44**(2), 489--514.
#'
#' @examples
#' pred <- PredPMS(twosines)
#'
#' @export

PredPMS <- function(data, modas, x.malla, h1 = 0.3, h2 = 0.5, dim.y = 1,
                    eps = 1e-8, conf.level = 0.95, k = 10, len = 200, malla){

  # comprobación de los argumentos suministrados:
  if(!methods::is(h1,"numeric") | h1<=0){
    stop("Bandwidth for covariates, `h1`, is not a positive (>0) numeric value.")
  }


  if(!methods::is(h2,"numeric")| h2<=0){
    stop("Bandwidth for response variable, `h2`, is not a positive (>0) numeric value.")

  }

  # calculamos la precisión
  p = -log(eps, base = 10)

  if(!methods::is(k,"numeric")){
    stop("The number of initial Y values per x.malla point must be an one-dimensional integer value.")
  }
  k = floor(abs(k))


  if (!methods::is(conf.level, "numeric") | length(as.vector(conf.level)) != 1L | conf.level > 1 | conf.level <0){
    stop("Confidence level must be a number between 0 and 1.")
  }


  if(!missing(malla)){
    if(!methods::is(malla, "MRpms_malla")) stop("`malla` is not a `malla` object. Use `mallador` function to compute one.")
  }

  # comprobar si recibimos una muestra válida
  if(missing(data)){
    stop("There is no sample data provided. Missing `data`")
  } else {

    if(methods::is(data,"data.frame")){ data <- as.matrix(data)}

    # comprobamos si de hecho es una matriz o array
    if((!methods::is(data,"array") & !methods::is(data,"matrix")) | typeof(data) != "double" ){
      stop("`data` argument is neither a numeric matrix nor an array.")

    }
  }


  # calculamos la dimensión de la covariable
  dim.x = ncol(data) - 1

  # Separamos la variable explicativa de la variable respuesta
  X = data[, 1:dim.x]
  Y = data[, dim.x+1:dim.y]

  # si no se da una malla, se especifica una
  if(missing(malla)){

    # de ser proveído un objeto modas, entonces construimos la malla con los puntos con los que esta fue calculada
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) {stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")}

      malla = mallador(data, x.malla = attr(modas,"x.malla"), k = k)
    } else{  # si no fue provisto un objeto `modas`, usamos los argumentos suministrados

      malla = mallador(data, k = k, len = len)
    }
  } else {

    # en caso de que tengamos ambos, comprobemos que están definidos sobre los mismos
    # puntos.
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")
      if(modas$x.malla != attr(malla,"x.malla")){
        stop("`x.malla` in `malla` object and `modas` object are not equal.
Please, provide a `malla` object built over the same `x.malla` as `modas`..")

      }
    }
  }

  # obtenemos la información necesaria de la malla
  x.malla <- attr(malla,"x.malla")
  k <- attr(malla, "k")
  len <- attr(malla, "len")
  n = length(Y)

  # en caso de que las modas no sean proporcionadas, se calculan.
  if(missing(modas)){
    modas <- PMSc(X = X, Y = Y, malla = malla, dim = dim.x,
                  h1 = h1, h2 = h2, p = p, eps = eps,
                  n = n, k = k, len = len)
    
    # la función PMSc es más rápida, pero no devuelve un objeto modas.
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
  }


  malla.aux = mallador(data, x.malla = X, k = k)
  # aux <- lapply(1:n, function(i) unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
  #                                                                  x = malla.aux[(i-1)*k+1,1:dim],
  #                                                                  ymalla = malla.aux[(i-1)*k+(1:k),dim+1],
  #                                                                  h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-2)))
  len.aux <- attr(malla.aux, "len")
  aux <- PMSc(X = X, Y = Y, malla = malla.aux, dim = dim.x,
                h1 = h1, h2 = h2, p = p, eps = eps,
                n = n, k = k, len = len.aux)


  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), conf.level)

  if(dim.x == 1L & dim.y == 1L){
    plot(X, Y)
    ni <- lapply(modas$modas, length)
    xg <- rep(x.malla, times = ni)
    yg <- unlist(modas$modas)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-epsh, yg[i]+epsh), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg, col = "red", pch = 19 )
  }
  salida <- list(modas = modas, epsh = epsh, conf.level = conf.level)
  class(salida) <- "MRpms_pred"
  return(salida)
}

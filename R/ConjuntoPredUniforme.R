#' Uniform prediction sets based on modal regression
#'
#' Computes a uniform prediction set for the response value.
#'
#'
#' @param muestra Matrix containing the sample.
#'
#' @param modas List containing the estimated set of conditional local
#' modes on each covariate point in `malla`. In case modas missing,
#' PMSc is called to compute them. It must be a MRpms_modas class object.
#'
#' @param x.malla Set of covariate points where the modes will be estimated.
#'  Not needed if `modas` is provided.
#'
#' @param h1 Smoothing parameter for covariates.
#'
#' @param h2 Smoothing parameter for response variable.
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
#' @return A list. First entry is a list containing the estimated set of
#' conditional local modes for each point in `x.malla`, the second element.
#' Third element is the radius of the prediction set.
#' Aditionally, the prediction set and the estimated modes are plotted.
#'
#' @references
#' Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L. (2016).
#' *Nonparametric modal regression*. The Annals of Statistics, **44**(2), 489--514.
#'
#' @examples
#' pred <- PredPMS(twosines)
#'
#' @export

PredPMS <- function(muestra, modas, x.malla, h1 = 0.3, h2 = 0.5,
                    eps = 1e-8, conf.level = 0.95, k = 10, len = 200, malla){

  # comprobación de los argumentos suministrados:
  if(!methods::is(h1,"numeric") | h1<=0){
    stop("Bandwidth for covariates, `h1`, is not a positive (>0) numeric value.")
  }


  if(!methods::is(h2,"numeric")| h2<=0){
    stop("Bandwidth for response variable, `h2`, is not a positive (>0) numeric value.")

  }

  # if(!methods::is(p,"numeric")){
  #   stop("Precision must be a one-dimensional numeric value.")
  # }
  # p = floor(abs(p))
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
  if(missing(muestra)){
    stop("There is no sample data provided. Missing `muestra` matrix")
  } else {


    # comprobamos si de hecho es una matriz o array
    if((!methods::is(muestra,"array") & !methods::is(muestra,"matrix")) | typeof(muestra) != "double" ){
      stop("`muestra` argument is not a numeric matrix nor an array.")

    }
  }


  # calculamos la dimensión de la covariable
  dim = ncol(muestra) - 1

  # comprobamos que el número de dimensiones sea correcto
  if (!methods::is(dim,"numeric")){
    stop("Number of dimensions is not correct. Please, check `muestra` has more
         than two columns.")

  }

  # Separamos la variable explicativa de la variable respuesta
  X = muestra[, 1:dim]
  Y = muestra[, dim+1]

  # si no se da una malla, se especifica una
  if(missing(malla)){

    # de ser proveído un objeto modas, entonces construimos la malla con los puntos con los que esta fue calculada
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) {stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")}
      # if (missing(k)){
      #   stop("Provide the desired number of Y values per x point on the `malla`, `k`.")
      #
      # }
      malla = mallador(muestra, x.malla = attr(modas,"x.malla"), k = k)
    } else{  # si no fue provisto un objeto `modas`, usamos los argumentos suministrados
      # if(missing(k) | missing(l)){
      #   stop("Not enough arguments to compute a `malla` object. Please, check `k` and `l` argument were provided.")
      #
      # }
      malla = mallador(muestra, k = k, len = len)
    }
  } else {
    # en caso de que tengamos ambos, comprobemos que están definidos sobre los mismos
    # puntos.
    if(!missing(modas)){
      if(!methods::is(modas, "MRpms_modas")) stop("`modas` is not an MRpms_modas object, please, use only an object result of an MRpms function package.")
      if(attr(modas,"x.malla") != attr(malla,"x.malla")){
        stop("`x.malla` in `malla` object and `modas` object are not equal.
Please, provide a `malla` object built over the same `x.malla` as `modas`..")

      }
    }
  }


  x.malla <- attr(malla,"x.malla")
  k <- attr(malla, "k")
  len <- attr(malla, "len")
  n = length(Y)

  if(missing(modas)){
        modas <- PMSc(X = X, Y = Y, malla = malla, dim = dim,
                      h1 = h1, h2 = h2, p = p, eps = eps,
                      n = n, k = k, len = len)
  }


  malla.aux = mallador(muestra, x.malla = X, k = k)
  # aux <- lapply(1:n, function(i) unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
  #                                                                  x = malla.aux[(i-1)*k+1,1:dim],
  #                                                                  ymalla = malla.aux[(i-1)*k+(1:k),dim+1],
  #                                                                  h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-2)))
  aux <- PMSc(X = X, Y = Y, malla = malla.aux, dim = dim,
                h1 = h1, h2 = h2, p = p, eps = eps,
                n = n, k = k, len = len)


  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), conf.level)

  if(dim == 1){
    plot(X, Y)
    ni <- lapply(modas, length)
    xg <- rep(x.malla, times = ni)
    yg <- unlist(modas)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-epsh, yg[i]+epsh), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg, col = "red", pch = 19 )
  }
  return(list(modas = modas, x.malla = x.malla, epsh = epsh))
  # return(epsh)
}

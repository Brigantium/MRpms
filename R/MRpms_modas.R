#' @title Plot method for `MRpms_modas` class with dimension `1`.
#' 
#' @description Plot method for `MRpms_modas` objects, which 
#' are computed from a `PMS` call or, a more bare-bones version, `PMSc` call. 
#' 
#' @param x An `MRpms_modas` object to plot.
#' 
#' @param data The data from which the modes were computed.
#' 
#' @param ... Other parameters
#' 
#' @returns A plot with both the `data` and the modes.
#' 
#' @seealso [PMS()]
#' 
#' 
#' @export
#' @method plot MRpms_modas
plot.MRpms_modas <- function(x, data, ...){
  if (all(x$dims == c(1L,1L))){
    plot(data)
    ni <- sapply(x$modas, length)
    xg <- rep(x$x.malla, times = ni)
    yg <- unlist(x$modas)
    graphics::points(xg, yg, ...)
  } else stop("This is not an univariate case in both response and covariable.")
}


#' @title Print an `MRpms_modas` object
#' 
#' @description A generic function to print an `MRpms_modas` object, which are computed
#' from a `PMS` cal or, a more bare-bones version, `PMSc` call.
#' 
#' @param x An `MRpms_modas` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return The modes in every covariable point from `x.malla`, which is an attribute of the `MRpms_modas`.
#' 
#' @export
#' @method print MRpms_modas

print.MRpms_modas <- function(x,...){
  y <- sapply(x$modas,as.numeric)
  print(y,...)
  # print(names(x))
  # NextMethod("print")
  # print.default(as.numeric(x))
  invisible(x)
}

#' @title Print useful information of an `MRpms_modas` object
#' 
#' @description A generic function to summarize useful information of an `MRpms_modas` object, which are computed
#' from a `PMS` call or, a more bare-bones version, `PMSc` call.
#' 
#' @param object An `MRpms_modas` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return The maximum, minimum and mean number of modes computed, with certain tolerance to avoid atipic values.
#' 
#' @export
#' @method summary MRpms_modas

summary.MRpms_modas <- function(object,...){
  ni <- sapply(object$modas,length)
  warning("To avoid atipic values, maximum and minimum number of modes are computed as 0.9 and 0.1 quantiles.")
  mean.modes <- mean(ni)
  cat(paste0("Mean of modes: ", mean.modes)); cat("\n")
  salida <- list(mean.modes = mean.modes)
  quantiles <- stats::quantile(ni,c(0.1,0.9))
  cat(paste0("Maximum number of modes in one point: ", quantiles[2])); cat("\n")
  cat(paste0("Minimum number of modes in one point: ", quantiles[1])); cat("\n")
  num.x <- length(object$x.malla)/object$dims[1]
  cat(paste0("Number of covariable points used: ", num.x)); cat("\n")
  salida <- append(salida, list(maximum.modes = quantiles[2], minimum.modes = quantiles[1], num.x = num.x))
  return(salida)
}
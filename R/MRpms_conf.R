#' @title Plot method for `MRpms_conf` class with dimension `1`.
#' 
#' @description Plot method for `MRpms_conf` objects, which 
#' are computed from a `ConfPMS` call and contains information to construct confidence intervals
#' for the computed modes. 
#' 
#' @param x An `MRpms_conf` object.
#' 
#' @param data The data from which the modes were computed.
#' 
#' @param ... Other parameters
#' 
#' @returns A plot with the `data`, in blank circles, the modes in red, the pointwise confidence intervals in blue
#' and the uniform confidence intervals in gray.
#' 
#' @seealso [ConfPMS()]
#' 
#' 
#' @export
#' @method plot MRpms_conf
plot.MRpms_conf <- function(x, data, ...){
  if (all(x$modas$dims == c(1L,1L))){
    plot(data)
    ni <- sapply(x$modas$modas, length)
    xg <- rep(x$modas$x.malla, times = ni)
    yg <- unlist(x$modas$modas)
    if (!is.null(x$delta)){
      delta <- x$delta
      sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-delta, yg[i]+delta), col ="lightgrey", lwd = 1.5))
      graphics::points(xg, yg - delta, col = "grey", pch = 19, cex = 0.4)
      graphics::points(xg, yg + delta, col = "grey", pch = 19, cex = 0.4)
    }
    if (!all(is.null(x$deltax))){
      radio <- rep(x$deltax, times = ni)
      sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-radio[i], yg[i]+radio[i]), col ="blue", lwd = 1.5))
    }
    
    graphics::points(xg, yg, col = "red", pch = 19)
  } else stop("This is not an univariate case in both response and covariable.")
}


#' @title Print an `MRpms_conf` object
#' 
#' @description A generic function to print an `MRpms_conf` object, which are computed
#' from a `ConfPMS` call.
#' 
#' @param x An `MRpms_conf` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return An index of contents.
#' 
#' @export
#' @method print MRpms_conf

print.MRpms_conf <- function(x,...){
  for(i in names(x)){
    cat(paste0("$",i));cat("\n")
  }
  # print(names(x),...)

  invisible(x)
}

#' @title Print useful information of an `MRpms_conf` object
#' 
#' @description A generic function to summarize useful information of an `MRpms_conf` object, which are computed
#' from a `ConfPMS` call.
#' 
#' @param object An `MRpms_conf` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return The summarize of the `MRpms_conf` object, plus estimated uniform error if computed and quantiles of pointwise errors if computed.
#' 
#' @export
#' @method summary MRpms_conf

summary.MRpms_conf <- function(object,...){
  sum.modas <- summary(object$modas)
  salida <- sum.modas
  # ni <- sapply(x$modas,length)
  # warning("To avoid atipic values, maximum and minimum number of modes are computed as 0.9 and 0.1 quantiles.")
  if(!is.null(object$delta)){
    cat(paste0("Uniform error: ", object$delta)); cat("\n")
    salida <- append(salida,list(unif.err = object$delta))
  }
  if(!is.null(object$deltax)){
    med.errx <- mean(object$deltax)
    cat(paste0("Mean of pointwise errors: ", med.errx)); cat("\n")
    cat("Relevant quantiles for pointwise errors\n")
    quantiles <- stats::quantile(object$deltax,c(0,0.25,0.5,0.75,1))
    print(quantiles)
    salida <- append(salida, list(med.errx = med.errx, quantiles.errx = quantiles))
  }
  invisible(salida)
}
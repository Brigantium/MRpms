#' @title Plot method for `MRpms_pred` class with dimension `1`.
#' 
#' @description Plot method for `MRpms_pred` objects, which 
#' are computed from a `PredPMS` call and contains information to construct prediction intervals
#' for the computed modes. 
#' 
#' @param x An `MRpms_pred` object.
#' 
#' @param data The data from which the modes were computed.
#' 
#' @param ... Other parameters
#' 
#' @returns A plot with the `data`, in blank circles and the uniform prediction intervals in gray.
#' 
#' @seealso [PredPMS()]
#' 
#' 
#' @export
#' @method plot MRpms_pred
plot.MRpms_pred <- function(x, data, ...){
  if (all(x$modas$dims == c(1L,1L))){
    plot(data)
    ni <- sapply(x$modas$modas, length)
    xg <- rep(x$modas$x.malla, times = ni)
    yg <- unlist(x$modas$modas)
    epsh <- x$epsh
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-epsh, yg[i]+epsh), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + epsh, col = "grey", pch = 19, cex = 0.4)
    
    graphics::points(xg, yg, col = "red", pch = 19)
  } else stop("This is not an univariate case in both response and covariable.")
}


#' @title Print an `MRpms_pred` object
#' 
#' @description A generic function to print an `MRpms_pred` object, which are computed
#' from a `PredPMS` call.
#' 
#' @param x An `MRpms_pred` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return The radius of the uniform prediction interval and the first modes.
#' 
#' @export
#' @method print MRpms_pred

print.MRpms_pred <- function(x,...){
  cat(paste0("Uniform radius for ",100*x$conf.level,"% of confidence: ", round(x$epsh,6))); cat("\n")
  # print(head(x$modas$modas))
  invisible(x)
}

#' @title Print useful information of an `MRpms_pred` object
#' 
#' @description A generic function to summarize useful information of an `MRpms_pred` object, which are computed
#' from a `PredPMS` call.
#' 
#' @param object An `MRpms_pred` object.
#' 
#' @param ... Other arguments. Currently, nothing more is avaliable.
#' 
#' @return The summarize of the `MRpms_modas` object, plus estimated radius for uniform prediction set and volume of said region.
#' 
#' @export
#' @method summary MRpms_pred

summary.MRpms_pred <- function(object,...){
  sum.modas <- summary(object$modas)
  salida <- sum.modas

  if (object$modas$dims[1] == 1L){
      lensopx <- max(object$modas$x.malla) - min(object$modas$x.malla)
    } else {
      lensopx <- sapply(1:object$modas$dims[1], function(i) max(object$modas$x.malla[,i])) - sapply(1:object$modas$dims[1], function(i) min(object$modas$x.malla[,i]))
    }
  
  ni <- sapply(object$modas$modas, length)
  intk <- mean(ni)*lensopx
  
  salida <- append(salida, list(vol.ph = intk*object$epsh))
  print(object)
  cat(paste0("Volume of the uniform prediction set: ",round(salida$vol.ph,5)))
  invisible(salida)
}
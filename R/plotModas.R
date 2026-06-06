#' @title Plot method for `MRpms_modas` class
#' 
#' @description 
#' 
#' 
#' @export
#' @method plot MRpms_modas
plot.MRpms_modas <- function(x, data, ...){
  if(!missing(data)){
    plot(data)
  } 
  # else if (attr(x,"data") != NULL) {
  #   plot(parse(attr(x,"data")))
  # }
  ni <- sapply(x, length)
  xg <- rep(attr(x, "x.malla"), times = ni)
  yg <- unlist(x)
  graphics::points(xg, yg, ...)
}

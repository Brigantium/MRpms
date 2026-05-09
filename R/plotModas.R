#' @exportS3Method

plot.MRpms_modas <- function(x, data, ...){
  if(!missing(data)){
    plot(data)
  }
  ni <- sapply(x, length)
  xg <- rep(attr(x, "x.malla"), times = ni)
  yg <- unlist(x)
  graphics::points(xg, yg, ...)
}

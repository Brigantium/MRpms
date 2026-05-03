# #' @exportS3Method
# plot.MRpms_modas <- function(x,...){
#   graphics::plot(X,Y)
#   ni <- lapply(x, length)
#   xg <- rep(attr(x,"x.malla"), times = ni)
#   yg <- unlist(x)
#   graphics::points(xg, yg, ...)
#   # NextMethod("plot")
# }
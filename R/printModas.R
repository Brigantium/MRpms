#' @exportS3Method

print.MRpms_modas <- function(x,...){
  y <- sapply(x,as.numeric)
  print(y,...)
  # NextMethod("print")
  # print.default(as.numeric(x))
  invisible(x)
}
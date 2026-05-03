#' @exportS3Method

print.MRpms_malla <- function(x, ...){
    y <- as.numeric(x)
  NextMethod("print")
  invisible(x)
}
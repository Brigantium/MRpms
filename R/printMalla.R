#' @export
#' @method print MRpms_malla

print.MRpms_malla <- function(x, ...){
  y <- as.numeric(x)
  print(y,...)
  # NextMethod("print")
  invisible(x)
}
#' @export 
## usethis namespace: start
#' @useDynLib MRpms, .registration = TRUE
## usethis namespace: end
  NULL


PMSc <- function(X,Y,malla = mallador(X,Y,dim = dim), dim = 1, h1, h2,p, eps, n = length(Y), k = attr(malla,"k"), l = attr(malla,"l")){

  .Call(PMS,as.numeric(X), as.numeric(Y), as.numeric(malla), as.integer(dim),as.numeric(h1),as.numeric(h2), as.integer(p), as.numeric(eps), as.integer(n),as.integer(k), as.integer(l))
}
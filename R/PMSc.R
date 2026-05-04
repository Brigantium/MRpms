  ## usethis namespace: start
  #' @useDynLib MRpms, .registration = TRUE
  ## usethis namespace: end
  NULL


PMSc <- function(X,Y,malla = mallador(cbind(X,Y)), dim = 1, h1, h2,p, eps, n = length(Y), k = attr(malla,"k"), len = attr(malla,"len")){

  .Call(PMSC,as.numeric(X), as.numeric(Y), as.numeric(malla), as.integer(dim),as.numeric(h1),as.numeric(h2), as.integer(p), as.numeric(eps), as.integer(n),as.integer(k), as.integer(len))
}

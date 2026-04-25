  ## usethis namespace: start
  #' @useDynLib MRpms, .registration = TRUE
  ## usethis namespace: end
  NULL

CVC <- function(X,Y, malla, dim, h1, h2, p, eps, n,k = attr(malla,"k")){

  .Call("CV", as.numeric(X), as.numeric(Y), as.numeric(malla), as.integer(dim),as.numeric(h1),as.numeric(h2),as.integer(p),as.numeric(eps),as.integer(n), as.integer(k))
}
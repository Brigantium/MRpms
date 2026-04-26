  ## usethis namespace: start
  #' @useDynLib MRpms, .registration = TRUE
  ## usethis namespace: end
  NULL

CVC <- function(muestra = cbind(X,Y),X = muestra[,1:dim],Y = muestra[,dim+1], 
                malla = mallador(X,Y, dim = dim), dim, h1, h2, p = floor(-log(eps, base = 10)),
                eps = 10^(-p), n = length(Y), k = attr(malla,"k")){

  .Call(CV, as.numeric(X), as.numeric(Y), as.numeric(malla), 
              as.integer(dim),as.numeric(h1),as.numeric(h2),
              as.integer(p),as.numeric(eps),as.integer(n), as.integer(k))
}
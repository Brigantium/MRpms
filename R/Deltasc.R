## usethis namespace: start
#' @useDynLib MRpms, .registration = TRUE
## usethis namespace: end
NULL

#  SEXP Xarg, SEXP Yarg, SEXP modasarg, SEXP nmodasarg, SEXP mallaarg, SEXP narg, SEXP karg, SEXP larg, SEXP h1arg, SEXP h2arg, SEXP parg, SEXP epsarg, SEXP dimarg, SEXP Barg, SEXP seed
Bdeltas <- function(X,Y,modas,ni,malla,n,k,l,h1,h2,p,eps,dim,B,seed){

  .Call(bootstrap_deltas,as.numeric(X), as.numeric(Y),as.numeric(unlist(modas)),as.integer(ni), as.numeric(malla),
          as.integer(n),as.integer(k), as.integer(l),as.numeric(h1), as.numeric(h2),as.integer(p),
        as.numeric(eps),as.integer(dim),as.integer(B),as.integer(seed))

}
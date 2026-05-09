Deltas <- function(X, Y, modas, malla, n = length(Y),
                   k = attr(malla,"k"), len = attr(malla,"len"),
                   h1  = 0.3, h2 = 0.5, p = floor(-log(eps, base = 10)),
                   eps = 10^(-p), dim){

  imb <- sample(1:n, size = n, replace = TRUE)
  Xb <- X[imb]
  Yb <- Y[imb]

  # modasb <- PMSc(Xb,Yb,malla = malla, dim = 1, h1 = h1, h2 = h2, p = p, eps = eps)
  modasb <- sapply(1:len, function(i) sort(unique(round(unlist(PMS1c(X = Xb, Y = Yb,
                                                                   x = malla[(i-1)*k+1,1:dim],
                                                                   ymalla = malla[(i-1)*k+(1:k),dim+1],
                                                                   h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-1))))
  Deltax <- sapply(1:len, function(i) .Haus(modas[[i]], modasb[[i]]))
  return(Deltax)
}

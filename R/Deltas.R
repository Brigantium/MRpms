.Deltas <- function(X, Y, modas, malla, n = length(Y),
                   k = attr(malla,"k"), l = attr(malla,"l"), 
                   h1  = 0.3, h2 = 0.5, p = floor(-log(eps, base = 10)), 
                   eps = 10^(-p), paso = .pasoU, dim, m){
  
  imb <- sample(1:n, size = n, replace = TRUE)
  Xb <- X[imb,]
  Yb <- Y[imb]


  modasb <- sapply(1:l, function(i) sort(unique(round(unlist(PMF1C(X = Xb, Y = Yb,
                                                                   x = malla[(i-1)*k+1,1:dim],
                                                                   ymalla = malla[(i-1)*k+(1:k),dim+1], 
                                                                   h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-1))))
  Deltax <- sapply(1:l, function(i) .Haus(modas[[i]], modasb[[i]]))
  return(Deltax)
}
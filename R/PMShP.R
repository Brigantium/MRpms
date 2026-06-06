PMS.hP <- function(data, dim = ncol(data)-1,
                   X = data[,1:dim], Y = data[,dim+1],
                   malla = mallador(data, x.malla = X),
                   hs.x = c(0.1,0.5,1), hs.y = c(0.1,0.5,1),
                   eps = 1e-8, p = floor(-log(eps, base = 10)),
                   conf.level = 0.95, ku = 10, lu = 100){

  
  Hs <- cbind(rep(hs.x, times = length(hs.y)),rep(hs.y, each = length(hs.x)))
  
  if (dim == 1L){
    lensopx <- max(X) - min(X)
  } else {
    lensopx <- sapply(1:dim, function(i) max(X[,i])) - sapply(1:dim, function(i) min(X[,i]))
  }
  n <- length(Y)
  

  medidashs <- apply(Hs, 1, function(h) hP(X = X, Y = Y, malla = malla,
                                            dim = dim, h1 = h[1], h2 = h[2],
                                            p = p, eps = eps, n = n,
                                            conf.level = conf.level, lensopx = lensopx))

  return(Hs[which.min(medidashs),])
}

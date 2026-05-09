PMS.hcv <- function(data, dim = ncol(data)-1,
                    X = data[,1:dim], Y = data[,dim+1],
                    malla = mallador(data, x.malla = X),
                    hs.x = c(0.1,0.5,1), hs.y = c(0.1,0.5,1),
                    eps = 1e-8, p = floor(-log(eps, base = 10))){

  # paso = if(dim==1) .pasoU else .pasoM
  Hs <- cbind(rep(hs.x, times = length(hs.y)),rep(hs.y, each = length(hs.x)))
  k <- attr(malla,"k")
  n <- length(Y)
  medidashs <- apply(Hs, 1, function(h) CVC(X = X, Y = Y, malla = malla,
                                            dim = dim, h1 = h[1], h2 = h[2],
                                            p = p, eps = eps,k = k, n = n))

  return(Hs[which.min(medidashs),])
}

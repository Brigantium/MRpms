PMS.hP <- function(muestra, dim = ncol(muestra)-1, 
                   X = muestra[,1:dim], Y = muestra[,dim+1], 
                   malla = mallador(X, Y, dim = dim, x.malla = X), 
                   hs.x = c(0.1,0.5,1), hs.y = c(0.1,0.5,1), 
                   eps = 1e-8, p = floor(-log(eps, base = 10)), 
                   nivel = 0.95, ku = 10, lu = 100){
  
  # paso = if(dim==1) .pasoU else .pasoM
  Hs <- cbind(rep(hs.x, times = length(hs.y)),rep(hs.y, each = length(hs.x)))
  
  n <- length(Y)
  mu <- nrow(mallau)
  
  medidashs <- apply(Hs, 1, function(h) hP(X = X, Y = Y, malla = malla, 
                                            dim = dim, h1 = h[1], h2 = h[2], 
                                            p = p, eps = eps, n = n, 
                                            nivel = nivel))
  
  return(Hs[which.min(medidashs),])
}
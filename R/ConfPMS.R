ConfPMS <- function(muestra, modas, dim = ncol(muestra) -1, 
                    X = muestra[,1:dim], Y = muestra[,dim+1], 
                    malla = mallador(X, Y, dim = dim), h1 = 0.3, h2 = 0.5, 
                    eps = 1e-8, p = floor(-log(eps, base = 10)),
                    nivel=0.95, B=500, type = 2,
                    semilla = 2026){

  
  xx <- attr(malla,"xx")
  k <- attr(malla, "k")
  l <- attr(malla, "l")
  m <- nrow(malla)
  
  if(missing(modas)){

    modas <- PMSc(X,Y,malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                  n = n, k = k, l = l)
  }
  plot(X, Y)
  ni <- sapply(modas, length)
  xg <- rep(xx, times = ni)
  yg <- unlist(modas)
  salida <- list(puntos=cbind(xg,yg))
  
  n <- length(Y)
  X <- matrix(X, ncol=dim)
  set.seed(semilla)
  
  Deltasbx <- replicate(B, .Deltas(X = X, Y = Y, modas = modas,
                                   malla = malla, n = n, k = k, l = l, 
                                   h1 = h1, h2 = h2, p = p, eps = eps, 
                                   dim = dim))
  
  if(type==1 | type==2){
    Deltasb <- apply(Deltasbx, 2, max)
    delta <- quantile(Deltasb, nivel)
    sapply(1:length(xg), function(i) lines(rep(xg[i],2), c(yg[i]-delta, yg[i]+delta), col ="lightgrey", lwd = 1.5))
    points(xg, yg - delta, col = "grey", pch = 19, cex = 0.4)
    points(xg, yg + delta, col = "grey", pch = 19, cex = 0.4)
    salida <- append(salida, list(delta = delta))
  }
  
  if(type==0 | type==2){
    deltasx <- apply(Deltasbx, 1, function(x) quantile(x, nivel))
    radio <- rep(deltasx, ni)
    # points(xg, yg - radio, col = "blue", pch = 19, cex = 0.4)
    # points(xg, yg + radio, col = "blue", pch = 19, cex= 0.4)
    sapply(1:length(xg), function(i) lines(rep(xg[i],2), c(yg[i]-radio[i], yg[i]+radio[i]), col ="blue", lwd = 1.5))
    salida <- append(salida, list(deltax = radio))
  }
  
  points(xg, yg, col = "red", pch = 19)
  
  return(salida)
}
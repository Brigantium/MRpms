#' @export

PredPMS <- function(muestra, modas, dim = ncol(muestra) -1, 
                                 X = muestra[,1:dim], Y = muestra[,dim+1], 
                                 malla = mallador(X, Y, dim = dim), 
                                 h1 = 0.4, h2 = 1, n = length(Y), 
                                 eps = 1e-8, p = -log(eps, base = 10), 
                                 nivel = 0.95){

  
  
  x.malla <- attr(malla,"x.malla")
  k <- attr(malla, "k")
  l <- attr(malla, "len")
  
  if(missing(modas)){
        modas <- PMSc(X,Y,malla = malla,h1 = h1,h2 = h2, p = p, eps = eps, dim = dim,
                  n = n, k = k, l = l)
  }

  malla.aux = mallador(X,Y,x.malla = X,dim = dim)
  aux <- sapply(1:n, function(i) unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
                                                                   x = malla.aux[(i-1)*k+1,1:dim],
                                                                   ymalla = malla.aux[(i-1)*k+(1:k),dim+1],
                                                                   h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-2)))
  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), nivel)
  
  if(dim == 1){
    plot(X, Y)
    ni <- sapply(modas, length)
    xg <- rep(x.malla, times = ni)
    yg <- unlist(modas)
    sapply(1:length(xg), function(i) graphics::lines(rep(xg[i],2), c(yg[i]-epsh, yg[i]+epsh), col ="lightgrey", lwd = 1.5))
    graphics::points(xg, yg - epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg + epsh, col = "grey", pch = 19, cex = 0.4)
    graphics::points(xg, yg, col = "red", pch = 19 )
  }
}
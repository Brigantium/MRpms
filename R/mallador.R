mallador <- function(X, Y, k = 10, l = 100, dim = 1, malla.x){

  if(dim > 1){
    
    if(missing(malla.x)){
      minx <- apply(X, 2, min)
      maxx <- apply(X, 2, max)
      xx <- apply(cbind(minx, maxx), 1, function(x) seq(x[1], x[2], len = l))
      malla.x <- sapply(1:dim, function(i) rep(rep(xx[,i], each = l^(dim-i)), times = l^(i-1)))
    }
    
    y.temp <- c(apply(malla.x, 1 ,function(x) Y[order(norma(x, X))[1:k]]))
    malla <- cbind(matrix(rep(malla.x, each = k), ncol = dim), y.temp)
  }else{
    
    if(missing(malla.x)){
      minx <- min(X)
      maxx <- max(X)
      malla.x <- seq(minx, maxx, len = l)
    } else {
      # malla.x <- unique(malla.x)
      l <- length(malla.x)
    }
    
    Xaux <- matrix(X,ncol=1)
    
    yg <- c(sapply(malla.x, function(x) Y[order(norma(x, Xaux))[1:k]]))
    xg <- rep(malla.x, each = k)
    malla <- cbind(xg,yg)
  }
  
  attr(malla,"l") <- l
  attr(malla,"k") <- k
  attr(malla, "xx") <- malla.x
  return(malla = malla)
}
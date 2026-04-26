PMSr <- function(muestra, malla = mallador(X,Y,dim = dim),
                dim, h1 = 0.5, h2 = 1, 
                p = floor(-log(eps,base = 10)), eps = 10^(-p),X = muestra[,1:dim],
                Y =muestra[,dim+1]){
  if(missing(dim)){
    message("Por favor, introdudza la dimension de X")
    return()
  }
  
  
  # paso = if(dim==1) .pasoU else .pasoM
  l = attr(malla, "len")
  k = attr(malla,"k")
  n = length(Y)
  
  modas <- sapply(1:l, function(i) sort(unique(round(unlist(PMS1c(X = X, Y = Y, x = malla[(i-1)*k +1,1:dim], 
                                                                              ymalla = malla[(i-1)*k+(1:k),dim+1], h1 = h1, 
                                                                              h2 = h2,eps = 1e-8,k = 10, n = n)),p-1))))
  attr(modas,"xx") <- attr(malla,"x.malla")
  return(modas)
}
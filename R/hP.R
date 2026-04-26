.hP<-function(X, Y, malla, dim = 1, h1, h2,  
              p = floor(-log(eps, base = 10)), eps = 10^(-p),
              n = length(Y), nivel = 0.95,
              mallau, ku, lu, pasoxxu, mu){
  
  k = attr(malla,"k")

  aux <- sapply(1:n, function(i) sort(unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
                                                                x = malla[(i-1)*k+1,1:dim],
                                                                ymalla = malla[(i-1)*k+(1:k),dim+1],
                                                                h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-1))))
  
  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), nivel)
  

  aux2 <- PMSc(X,Y,malla = mallau,h1 = h1,h2 = h2, p = p, dim = dim,
                l = lu, k = ku)
  
  ni <- sapply(aux2, length)
  intk <- sum(ni)*pasoxxu

  
  return(epsh*intk)
}
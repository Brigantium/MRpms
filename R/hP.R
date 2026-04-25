.hP<-function(X, Y, malla, dim = 1, h1, h2,  
              p = floor(-log(eps, base = 10)), eps = 10^(-p),
              n = length(Y), paso = .pasoU, nivel = 0.95,
              mallau, ku, lu, pasoxxu, mu){
  
  k = attr(malla,"k")
  # aux <- .PMS2CV(X = X, Y = Y, malla = malla, h1 = h1, h2 = h2, p = p, eps=eps,
  #                paso = paso, dim = dim)
  aux <- sapply(1:n, function(i) sort(unique(round(unlist(PMF1C(X = X[-i], Y = Y[-i],
                                                                x = malla[(i-1)*k+1,1:dim],
                                                                ymalla = malla[(i-1)*k+(1:k),dim+1],
                                                                h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-1))))
  epsh <- quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), nivel)
  
  # aux2 <- .PMS2(X, Y, malla = mallau, h1 = h1, h2 = h2, p = p, eps = eps,
  #               l = lu, k = ku, paso = paso, dim = dim, m = mu)
  aux2 <- PMSc(X,Y,malla = mallau,h1 = h1,h2 = h2, p = p, dim = dim,
                paso = paso, l = lu, k = ku)
  
  ni <- sapply(aux2, length)
  intk <- sum(ni)*pasoxxu

  
  return(epsh*intk)
}
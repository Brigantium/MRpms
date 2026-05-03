hP <-function(X, Y, malla, dim = 1, h1, h2,
              p = floor(-log(eps, base = 10)), eps = 10^(-p),
              n = length(Y), nivel = 0.95, k = attr(malla, "k"),
              lensopx = max(X) - min(X)){

  aux <- lapply(1:n, function(i) sort(unique(round(unlist(PMS1c(X = X[-i], Y = Y[-i],
                                                                x = malla[(i-1)*k+1,1:dim],
                                                                ymalla = malla[(i-1)*k+(1:k),dim+1],
                                                                h1 = h1, h2 = h2,eps = eps, k = k, n = n)),p-2))))
  epsh <- stats::quantile(sapply(1:n, function(i) min(abs(Y[i] - aux[[i]]))), nivel)

  ni <- sapply(aux, length)
  intk <- mean(ni)*lensopx

  return(epsh*intk)
}
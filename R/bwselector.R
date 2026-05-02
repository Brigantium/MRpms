bwselector <- function(muestra, dim = ncol(muestra)-1,
                       X = muestra[,1:dim], Y = muestra[,dim+1],
                       malla = mallador(X, Y, dim = dim, x.malla = X),
                       eps = 1e-8, p = floor(-log(eps, base = 10)),
                       nivel = 0.95, metodo="CV", n = length(Y), 
                       eps2 = 1e-3){

  if(metodo == "P"){
    k <-attr(malla, "k")
    lensopx <- max(X) - min(X)
    medida <- function(h){hP(X = X, Y = Y, malla = malla, dim = dim,
                                  h1 = h[1], h2 = h[2], eps = eps, p = p,
                                  n = n, nivel = nivel, k = k, 
                                  lensopx = lensopx)}
    
  }else{
    medida <- function(h){CVC(X = X, Y = Y, malla = malla, dim = dim,
                              h1 = h[1], h2 = h[2], p = p, eps = eps, n = n)}
  }

  nsimplex <- dim + 1
  h0 <- c(max(X)-min(X),max(Y)-min(Y))/10
  as <- min(h0)
  qs <- as/(nsimplex)/sqrt(2)*(sqrt(nsimplex+1)-1)
  hs <- t(rbind(h0, matrix(rep(h0 + rep(qs,nsimplex),2),byrow=TRUE,nrow=nsimplex) + diag(as/sqrt(2),nsimplex)))
  
  # plot(t(hs), pch = 19, col = "red", xlim=c(0,h0[1]+3*as),ylim=c(0,h0[2]+3*as)) ###
  # lines(t(cbind(hs,hs[,1])), col = "red") ###
  
  auxX <- as.matrix(dist(X)) + max(dist(X))*diag(n)
  auxY <- as.matrix(dist(Y)) + max(dist(Y))*diag(n)
  hminX <- max(apply(auxX,1,min))
  hminY <- max(apply(auxY,1,min))
  htol <- c(hminX, hminY)/2

  lhs <- log(hs-htol)
  
  # for(i in 1:ncol(hs)){cat("Probando h=",hs[,i],"\n")} ###
  medidashs <- apply(hs, 2, medida)
  
  lhb0 <- lhs[,which.min(medidashs)]

  it <- 0
  j <- 0
  
  fbueno <- medidashs[which.min(medidashs)]
  while(sd(medidashs)/fbueno > eps2){

    orden <- order(medidashs, decreasing = TRUE)
    bueno <- orden[nsimplex+1]
    feo <- orden[2]
    malo <- orden[1]

    xmalo <- lhs[,malo]

    fbueno <- medidashs[bueno]
    ffeo <- medidashs[feo]
    fmalo <- medidashs[malo]

    x0 <- rowMeans(lhs[,-malo])

    xr <- 2*x0 - xmalo
    # cat("Probando h=",exp(xr)+htol,"\n") ###
    fxr <- medida(exp(xr)+htol)

    if(fxr >= ffeo){ # Contraction
      if(fxr < fmalo){
        xc <- (x0+xr)/2
        # cat("Probando h=",exp(xc)+htol,"\n") ###
        fxc <- medida(exp(xc)+htol)
        if(fxc < fxr){
          lhs[,malo] <- xc
          medidashs[malo] <- fxc
        }else{ # Shrink
          xbueno <- lhs[,bueno]
          lhs[,-bueno] <- apply(lhs[,-bueno], 2, function(x) (x + xbueno)/2)
          # for(i in (1:ncol(lhs))[-bueno]){
          #   cat("Probando h=",exp(lhs[,i])+htol,"\n") ###
          # }
          medidashs[-bueno] <- apply(exp(lhs[,-bueno])+htol, 2, function(h) medida(h))
        }
      }else{
        xc <- (x0 + xmalo)/2
        # cat("Probando h=",exp(xc)+htol,"\n") ###
        fxc <- medida(exp(xc)+htol)
        if(fxc < fmalo){
          lhs[,malo] <- xc
          medidashs[malo] <- fxc
        }else{ # Shrink
          xbueno <- lhs[,bueno]
          lhs[,-bueno] <- apply(lhs[,-bueno], 2, function(x) (x + xbueno)/2)
          # for(i in (1:ncol(lhs))[-bueno]){
          #   cat("Probando h=",exp(lhs[,i])+htol,"\n") ###
          # }
          medidashs[-bueno] <- apply(exp(lhs[,-bueno])+htol, 2, function(h) medida(h))
        }
      }
    }else{
      if(fxr >= fbueno){ # Reflection
        lhs[,malo] <- xr
        medidashs[malo] <- fxr
      }else{ # Expansion
        xe <- 2*xr - x0
        # cat("Probando h=",exp(xe)+htol,"\n") ###
        fxe <- medida(exp(xe)+htol)
        if(fxe < fxr){
          lhs[,malo] <- xe
          medidashs[malo] <- fxe
        }else{
          lhs[,malo] <- xr
          medidashs[malo] <- fxr
        }
      }
    }
    
    # plot(t(exp(lhs)+htol), pch = 19, col = "red", xlim=c(0,h0[1]+3*as),ylim=c(0,h0[2]+3*as)) ###
    # lines(t(exp(cbind(lhs,lhs[,1]))+htol), col = "red") ###

    it=it+1
    if(it == 100){break}

    lhb <- lhs[,which.min(medidashs)]
    if(sqrt(sum(((exp(lhb)-exp(lhb0))/exp(lhb))^2)) < eps2){j <- j+1}else{j <- 0}
    if(j == 10){break}
    lhb0 <- lhb
    
    # cat("Desv. Típica=",sd(medidashs),", it=",it,", j=",j,"\n") ###
    
  }

  return(exp(lhs[,which.min(medidashs)])+htol)
}
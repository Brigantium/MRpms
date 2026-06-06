#' @title Bandwith selector
#'
#' @description Performs an heuristic search based on the Nedler-Mead algorithm to find
#' appropiate values for the bandwith parameters `h1` and `h2`.
#' (Warning: this is a computationally requiring method which
#' might take some time to finish).
#'
#' @param data Matrix or data frame containing the sample. Currently, the response values are needed to be in
#' the last column.
#'
#' @param malla Matrix containing the initial points used to compute the modes
#' using the Partial Mean-Shift algorithm. If not provided,
#' `mallador`function is called. It must be a MRpms_malla class object.
#'
#' @param method `CV` (default) for Zhao and Huang selector or `P` for
#' Chen et al. selector based on the size of prediction sets.
#'
#' @param conf.level Confidence level for prediction sets (if `metodo`=`P`).
#'
#' @param eps Convergence tolerance for Partial Mean-Shift algorithm.
#'
#' @param eps2 Convergence tolerance for heuristic algorithm.
#' 
#' @param maxiter Maximum number of iterations for heuristic algorithm.
#' 
#' @param maxiter.tol Maximum number of consecutive iterations without improvement in heuristic algorithm before termination.
#'
#' @return
#' A vector with two elements, `h1` and `h2`, first one been for covariable and second one for response.
#'
#' @references
#' Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L. (2016).
#' *Nonparametric modal regression*. The Annals of Statistics, **44**(2), 489--514.
#'
#' Zhou, H. and Huang, X. (2019).
#' *Bandwidth selection for nonparametric modal regression*.
#' Communication in Statistics - Simulation and Computation,
#' **48**(4), 968--984.
#' 
#' Luersen, M., Le Riche, R. and Guyon, F. (2004). 
#' _A constrained, globalized, and bounded Nelder–Mead method for engineering optimization._
#' Structural and Multidisciplinary Optimization, **27**(1). 43--54. 
#' 
#' 
#'
#' @examples
#' h <- bwselector(twosines)
#' modas <- PMS(twosines, h1 = h[1], h2 = h[2])
#'
#' plot(twosines)
#' plot(modas, pch = 19, col = "red")
#'
#' @export

bwselector <- function(data, malla = mallador(data, x.malla = X),
                       eps = 1e-8, conf.level = 0.95, method="CV", eps2 = 1e-3,
                       maxiter = 100, maxiter.tol = 10){

  # calculamos la dimensión de la covariable
  dim <- ncol(data) - 1

  # separamos la covariable de la variable respuesta
  X <- data[,1:dim]
  Y <- data[,dim+1]

  # calculamos el tamaño muestral
  n <- length(Y)

  # calculamos la precisión 
  p <- floor(-log(eps, base = 10))

  # elegimos el método:
  # --------------------
  # si method = "P" entonces empleamos el volumen de los conjuntos de predicción
  # si method = "CV" entonces se emplea la extensión de la cross-validation 
  #
  if(method == "P"){
    k <-attr(malla, "k")

    # cálculo del volumen del soporte de la covariable
    if (dim == 1L){
      lensopx <- max(X) - min(X)
    } else {
      lensopx <- sapply(1:dim, function(i) max(X[,i])) - sapply(1:dim, function(i) min(X[,i]))
    }

    # definimos entonces la medida como el voluen de los conjuntos de predicción 
    medida <- function(h){hP(X = X, Y = Y, malla = malla, dim = dim,
                                  h1 = h[1], h2 = h[2], eps = eps, p = p,
                                  n = n, conf.level = conf.level, k = k,
                                  lensopx = lensopx)}

  }else{

    # en el caso de la cross-validation empleamos la función sugerida por Zhou y Huang.
    medida <- function(h){CVC(X = X, Y = Y, malla = malla, dim = dim,
                              h1 = h[1], h2 = h[2], p = p, eps = eps, n = n)}
  }

  # definimos el número de vértices, es decir el número de h buscados más uno
  nsimplex <- dim + 1
  h0 <- c(max(X)-min(X),max(Y)-min(Y))/10 # h inicial (va a ser necesario cambiarlo)
  as <- min(h0) # volumen inicial de simplex
  qs <- as/(nsimplex)/sqrt(2)*(sqrt(nsimplex+1)-1) # ver paper A constraines, globalized and bound Nelder-Mead
  hs <- t(rbind(h0, matrix(rep(h0 + rep(qs,nsimplex),2),byrow=TRUE,nrow=nsimplex) + diag(as/sqrt(2),nsimplex))) # los otros h para construir el simplex inicial

  # plot(t(hs), pch = 19, col = "red", xlim=c(0,h0[1]+3*as),ylim=c(0,h0[2]+3*as)) ###
  # lines(t(cbind(hs,hs[,1])), col = "red") ###

  # invención de ángel para que no podamos tener ceros en el denominador del mean shift
  auxX <- as.matrix(stats::dist(X)) + max(stats::dist(X))*diag(n)
  auxY <- as.matrix(stats::dist(Y)) + max(stats::dist(Y))*diag(n)
  hminX <- max(apply(auxX,1,min))
  hminY <- max(apply(auxY,1,min))
  htol <- c(hminX, hminY)/2
  #------------

  # transformación por ser el espacio de las h acotado
  lhs <- log(hs-htol)

  # for(i in 1:ncol(hs)){cat("Probando h=",hs[,i],"\n")} ###

  # evaluación de los h en la medida de error
  medidashs <- apply(hs, 2, medida)


  # tomamos el mejor
  bueno <- which.min(medidashs)
  lhb0 <- lhs[,bueno]


  it <- 0 # contador de iteraciones
  j <- 0 # número de iteraciones sin mejorar

  fbueno <- medidashs[bueno]  # evaluación de la medida en el mejor

  while(stats::sd(medidashs)/fbueno > eps2){ # condición de parada: que la desviación típica reescalada por el mejor de los puntos sea menor a un cierto épsilon pes2

    # buscamos el mejor, el peor y el segundo peor (feo)
    orden <- order(medidashs, decreasing = TRUE)
    bueno <- orden[nsimplex+1]
    feo <- orden[2]
    malo <- orden[1]

    # el punto en si asociado al malo
    xmalo <- lhs[,malo]

    # los valores de dichos puntos
    fbueno <- medidashs[bueno]
    ffeo <- medidashs[feo]
    fmalo <- medidashs[malo]

    # calculamos el punto medio de todos los puntos menos el peor
    x0 <- rowMeans(lhs[,-malo])

    xr <- 2*x0 - xmalo
    # cat("Probando h=",exp(xr)+htol,"\n") ###
    fxr <- medida(exp(xr)+htol)

    if(fxr >= ffeo){             # Contraction

      if(fxr < fmalo){
        xc <- (x0+xr)/2

        # cat("Probando h=",exp(xc)+htol,"\n") ###

        fxc <- medida(exp(xc)+htol)
        if(fxc < fxr){
          lhs[,malo] <- xc
          medidashs[malo] <- fxc

        }else{                   # Shrink
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

        }else{                 # Shrink
          xbueno <- lhs[,bueno]
          lhs[,-bueno] <- apply(lhs[,-bueno], 2, function(x) (x + xbueno)/2)

          # for(i in (1:ncol(lhs))[-bueno]){
          #   cat("Probando h=",exp(lhs[,i])+htol,"\n") ###
          # }

          medidashs[-bueno] <- apply(exp(lhs[,-bueno])+htol, 2, function(h) medida(h))
        }
        
      }
    }else{
      if(fxr >= fbueno){            # Reflection

        lhs[,malo] <- xr
        medidashs[malo] <- fxr

      }else{                        # Expansion

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

    it=it+1 # aumentamos el número de iteraciones
    if(it == maxiter) break # en caso de llegar al número máximo de iteraciones terminamos

    lhb <- lhs[,which.min(medidashs)]

    if(sqrt(sum(((exp(lhb)-exp(lhb0))/exp(lhb))^2)) < eps2) j <- j+1 else j <- 0  # si mejoramos, reiniciamos el contador

    if(j == maxiter.tol) break # en caso de llevar muchas iteraciones sin mejorar, terminamos 
    lhb0 <- lhb

    # cat("Desv. Típica=",sd(medidashs),", it=",it,", j=",j,"\n") ###

  }

  return(exp(lhs[,which.min(medidashs)])+htol)
}

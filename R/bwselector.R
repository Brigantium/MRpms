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
#' @param dim.y Response dimension.
#'
#' @param method `CV` (by default) for Zhou and Huang selector or `P` for
#' Chen et al. selector based on the size of prediction sets.
#'
#' @param conf.level Confidence level for prediction sets (if `method = P`).
#'
#' @param eps Convergence tolerance for Partial Mean-Shift algorithm.
#'
#' @param eps2 Convergence tolerance for heuristic algorithm.
#' 
#' @param maxiter Maximum number of iterations for heuristic algorithm.
#' 
#' @param maxiter.tol Maximum number of consecutive iterations without improvement in heuristic algorithm before termination.
#'
#' @param alpha Parameter of reflection in the Nelder-Mead algorithm.
#' 
#' @param gamma Parameter of expansion in the Nelder-Mead algorithm.
#' 
#' @param rho Parameter of contraction in the Nelder-Mead algorithm.
#' 
#' @param sigma Parameter of shrink in the Nelder-Mead algorithm.
#' 
#' @param initial.h Candidates to construct the initial simplex in the Nelder-Mead algorithm. Currently, there are needed
#' only three diferent candidates saves by rows.
#' 
#' 
#' @return
#' A vector with two elements, `h1` and `h2`, first one been for covariable and second one for response.
#'
#' @details
#' This function allows a Nelder-Mead based search of smooth parameter. There is two diferent measurments: 
#' 
#' * If `method = CV` then the algorithm uses the cross-validation like formula suggested by Zhou and Huang as degree of 
#' fitting,
#' \deqn{f(\mathbf h) = \displaystyle \frac{1}{n}\sum_{i=1}^n d^2\left(\widehat{M}_{n,\mathbf{h},-i}\left(X_i\right),Y_i\right)\widehat{N}_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right),}
#' where \eqn{\widehat M_{n,\mathbf{h},-i}(X_i)} is the estimated modes in \eqn{X_i} without consider \eqn{X_i} and \eqn{\widehat{N}_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right)}
#' is the number of modes estimated in \eqn{X_i}.
#' 
#' * If `method = P` then uses the volume of prediction sets instead as measurement of fitting,
#' \deqn{f(\mathbf h) = \displaystyle Vol\left(\widehat{\mathcal{P}}_{1-\alpha,\mathbf{h}}\right)=\widehat{\varepsilon}_{1-\alpha,\mathbf{h}}\int_{x\in D}\widehat{N}_{\mathbf{h}}(x)\, dx,}
#' where \eqn{\mathcal {\widehat P}_{1-\alpha, \mathbf {h}}} is the estimation of the prediction set, i.e. a set for which \eqn{\mathbb{P}(Y\in \mathcal P \geq 1-\alpha)};
#' \eqn{\hat \varepsilon_{1-\alpha,\mathbf h}} is the estimation of the \eqn{1-\alpha} quantile of the random
#' variable \eqn{d(Y,M(X))}, computed as the sample quantile of the sample \eqn{\{d(Y_i,\widehat{M}_{n,\mathbf h}(X_i))\,\colon i \in \{1,\dots,n\}\}};
#' and \eqn{\widehat{N}_{\mathbf h}(x)} is the estimated number of modes condicional to \eqn{X = x}.
#' 
#' ### **Nelder-Mead algorithm:** 
#' 
#' Selected the ranking method, the algorithm preforms a Nelder-Mead heuristic strategy to minimize \eqn{f}. This proccess 
#' starts with 3 different candidates to \eqn{\mathbf h}: \eqn{\mathbf h_1}, \eqn{\mathbf h_2} and \eqn{\mathbf h_3}. 
#' In each one computes \eqn{f(\mathbf h_i)}. Now, the idea is the same in each iteration: 
#' 
#' 1. Firstly, the \eqn{f(\mathbf h_i)} values are ordered, without lost of generality: \eqn{f(\mathbf h_1) \leq f(\mathbf h_2) \leq f(\mathbf h_3)}.
#' 
#' 2. the centroid of \eqn{\mathbf h_1} and \eqn{\mathbf h_2} is computed (\eqn{\mathbf h_0}), i.e., the mean point in coordinate terms —plus certain tolerance to avoid zeros—. 
#' 
#'     NOTE: _the original Nelder-Mead algorithm allows
#'         only unrestricted spaces, in other words the entirety of \eqn{\mathbb R^n}. To bypass this, we previously transform through 
#'         a natural logarithm function the espace where the smooth parameters lives, \eqn{(\mathbb{R}^{+})^n}, to \eqn{\mathbb R^n}. So, we 
#'         will work always with \eqn{\log \mathbf h} except when evaluate \eqn{f}._
#'
#' 3. **Reflection**: we compute a new candidate: \eqn{\mathbf h_r = \mathbf h_0 + \alpha(\mathbf h_0 - \mathbf h_3)}, with \eqn{\alpha > 0}, and compute\eqn{f(\mathbf h_r)}:
#' 
#'     * if \eqn{f(\mathbf h_1) < f(\mathbf h_r) < f(\mathbf h_2)} (not enough improvement): replace \eqn{\mathbf h_3} with \eqn{\mathbf h_r} and go to (1.).
#' 
#'     * if \eqn{f(\mathbf h_r) < f(\mathbf h_1)} (improvement): go to (4.).
#' 
#'     * if \eqn{f(\mathbf h_2) <f(\mathbf h_r)} (no improvement): go to (5.).
#' 
#'     NOTE: _we use \eqn{\alpha = 1}_.
#'  
#' 4. **Expansion**: we compute a new candidate: \eqn{\mathbf h_e = \mathbf h_0 + \gamma(\mathbf h_r - \mathbf h_0)}, with \eqn{\gamma  >1}, then
#' 
#'     * if \eqn{f(\mathbf h_e) < f(\mathbf h_r)} then replace \eqn{\mathbf h_3} with \eqn{\mathbf h_e} and go to (1.).
#' 
#'     * if \eqn{f(\mathbf h_r) < f(\mathbf h_e)} then replace \eqn{\mathbf h_3} with \eqn{\mathbf h_r} and go to (1.).
#' 
#'     NOTE: _we use \eqn{\gamma = 2}_.
#' 
#' 5. **Contraction**: two situations:
#' 
#'     * if \eqn{f(\mathbf h_r) < f(\mathbf h_3)} then compute a new candidate: \eqn{\mathbf h_c = \mathbf h_0 + \rho (\mathbf h_r - \mathbf h_0)}, with \eqn{0<\rho \leq 0.5}.
#' 
#'     * if \eqn{f(\mathbf h_3) < f(\mathbf h_r)} then compute a new condadte: \eqn{\mathbf h_c = \mathbf h_0 + \rho (\mathbf h_3 - \mathbf h_0)}, with \eqn{0 < \rho \leq 0.5}.
#' 
#'     In case of improvement, i.e. \eqn{f(\mathbf h_c) < \min\{f(\mathbf h_3),f(\mathbf h_r)\}},  replace \eqn{\mathbf h_3} with \eqn{\mathbf h_c}. In other case, go to (6.).
#'  
#'     NOTE: _we use \eqn{\rho = 1/2}_.
#' 
#' 6. **Shrink**: Fixed \eqn{\mathbf h_1}, replace the others with \eqn{\mathbf h_i = \mathbf h_1 + \sigma (\mathbf h_i - \mathbf h_1)}. Returns to (1.).
#' 
#'     NOTE: _we take \eqn{\sigma = 1/2}_.
#' 
#'  **Exit criteria:** The heuristic strategy ends when anyone of this statements is true:
#' 
#'  * The algorithm performs `maxiter` iterations or there is no improvement in `maxiter.tol` consecutive iterations.
#' 
#'  * The standard deviation of \eqn{f(\mathbf h_1),\ f(\mathbf h_2)} and \eqn{f(\mathbf h_3)}, reescalated by the value of \eqn{f(\mathbf h_1)} is less than a certain tolerance `eps`.
#' 
#' Further details could be reviewed in Luersen, Le Riche and Guyon (2004).
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
#' plot(modas,twosines, pch = 19, col = "red")
#'
#' @export


bwselector <- function(data, malla = mallador(data, x.malla = X), dim.y = 1,
                       eps = 1e-8, conf.level = 0.95, method="CV", eps2 = 1e-3,
                       maxiter = 100, maxiter.tol = 10, alpha = 1, gamma = 2, rho = 1/2, sigma = 1/2, initial.hs){

  # calculamos la dimensión de la covariable
  dim.x <- ncol(data) - dim.y

  # separamos la covariable de la variable respuesta
  X <- data[,1:dim.x]
  Y <- data[,dim.x+dim.y]

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
    if (dim.x == 1L){
      lensopx <- max(X) - min(X)
    } else {
      lensopx <- sapply(1:dim.x, function(i) max(X[,i])) - sapply(1:dim.x, function(i) min(X[,i]))
    }

    # definimos entonces la medida como el voluen de los conjuntos de predicción 
    medida <- function(h){hP(X = X, Y = Y, malla = malla, dim = dim.x,
                                  h1 = h[1], h2 = h[2], eps = eps, p = p,
                                  n = n, conf.level = conf.level, k = k,
                                  lensopx = lensopx)}

  }else{

    # en el caso de la cross-validation empleamos la función sugerida por Zhou y Huang.
    medida <- function(h){CVC(X = X, Y = Y, malla = malla, dim = dim.x,
                              h1 = h[1], h2 = h[2], p = p, eps = eps, n = n)}
  }

  
  nsimplex <- dim.x + dim.y

  if(missing(initial.hs)){
    h0 <- c(max(X)-min(X),max(Y)-min(Y))/10 # h inicial (va a ser necesario cambiarlo)
    as <- min(h0) # volumen inicial de simplex
    qs <- as/(nsimplex)/sqrt(2)*(sqrt(nsimplex+1)-1) # ver paper A constraines, globalized and bound Nelder-Mead
    hs <- t(rbind(h0, matrix(rep(h0 + rep(qs,nsimplex),2),byrow=TRUE,nrow=nsimplex) + diag(as/sqrt(2),nsimplex))) # los otros h para construir el simplex inicial
  } else {

    if (is.null(dim(initial.hs))) stop("Please, provided an array with one candidate per row.")
    
    n.aux <- nrow(initial.hs)

    if(n.aux < nsimplex+1) stop(paste0("Please, provided enough candidates, for this problem: ",nsimplex+1))
    
    alldiferent <- TRUE
    # chek if there is any candidate equal to other
    for(i in 1:(n.aux-1)){
      candidate.aux <- initial.hs[i,]
      for(j in (i+1):n.aux){
        if(all(abs(candidate.aux - initial.hs[j,]) < eps2)){
          alldiferent <- FALSE
          stop("Please, provided different candidates.")
        }
      }
    }

    if(alldiferent){
       hs <- t(initial.hs)
    } else{
      stop("Please, provided a valid initial set of candidates.")
    }
  }


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



    # Step 1: evaluation
    #####################################################
    # we search for 
    # - the best candidate (bueno), 
    # - the worst (malo)
    # - and the second worst (feo). 

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

    #- end: evaluaton -#################################

    # --------------------------------------------------

    # Step 2: centroid
    ####################################################
    # The centroid is computed using all the candidates save
    # the worst (malo).

    # calculamos el punto medio de todos los puntos menos el peor
    x0 <- rowMeans(lhs[,-malo])

    #- end: centroid -##################################

    # --------------------------------------------------

    # Step 3: reflection
    ####################################################
    # A new candidate is computed, h_r, in the opposite
    # direction to h_3 - h_0. 
    #
    #         h_r <··········· h_0 --------------- h_3
    # Proportion:    alpha

    # xr <- 2*x0 - xmalo  # versión optimizada si usamos los parámetros usuales de  # reflection

    xr <- x0 + alpha*(x0 - xmalo)              # reflection

    # cat("Probando h=",exp(xr)+htol,"\n") ###

    fxr <- medida(exp(xr)+htol)    # evaluamos en el nuevo punto

    #- end: reflection #################################


    # now we check: f(h_r) > f(h_feo)
    # - true: we proceded with Step 5: Contraction
    # - false: we proceded with Step (4.) or (1.)
    if(fxr >= ffeo){              

      # --------------------------------------------------------------------

      # check if the new candidate is better than the worst, f(h_r) < f(h_malo),
      # - true: we proceded with Step 5: Contraction from the reflection candidate, h_r.
      # - false: we proceded with Step 5: Contraction from the worst candidate, h_3.
      if(fxr < fmalo){

        # Step 5: Contraction (case f(h_r) < f(h_malo))
        ####################################################################
        # If f(h_r) < f(h_malo) then we contract the new candidate through 
        # the segment [h_0,h_r],
        #                
        #          h_0 ------- h_c <······················· h_r
        # Proportion:    rho               1 - rho

        # xc <- (x0+xr)/2  # optimized version when the usal parameters are used of   # Contraction

        xc <- x0 + rho*(xr-x0)    # contraction

        # cat("Probando h=",exp(xc)+htol,"\n") ###

        fxc <- medida(exp(xc)+htol)   # f(h_c)


        if(fxc < fxr){  

          # in case of not worsening, replace h_malo with h_c
          lhs[,malo] <- xc
          medidashs[malo] <- fxc

        }else{              # en caso contrario, pasamos al Paso 6

          # Paso 6: Shrink 
          #######################################################
          # we shrink the simplex fixing the best candidate.
          # 

          xbueno <- lhs[,bueno]
          # lhs[,-bueno] <- apply(lhs[,-bueno], 2, 
          # function(x) (x + xbueno)/2)    # optimized version
                                           #  when the usual parameters
                                           #  are used                  # shrink

          lhs[,-bueno] <- apply(lhs[,-bueno], 2, 
            function(x) xbueno + sigma*(x - xbueno)/2)            # shrink
          
          # for(i in (1:ncol(lhs))[-bueno]){
          #   cat("Probando h=",exp(lhs[,i])+htol,"\n") ###
          # }

          medidashs[-bueno] <- apply(exp(lhs[,-bueno])+htol, 2, function(h) medida(h))

          #- end: Shrink -#########################################

        }

        #- end: contraction -################################################


      }else{

        # Step 5: Contraction (case f(h_r) > f(h_malo))
        ###################################################
        # If f(h_r) > f(h_malo) then we compute through the
        # segment [h_0,h_malo],
        #
        #         h_0 ------- h_c <························ h_malo
        # proportion:   rho               1 - rho

        # xc <- (x0 + xmalo)/2    # optimized version when the usual parameters are used     # contraction

        xc <- x0 + rho*(xmalo - x0)  # contraction

        # cat("Probando h=",exp(xc)+htol,"\n") ###
        
        fxc <- medida(exp(xc)+htol)

        if(fxc < fmalo){

          # in case of not worsening, replace h_malo with h_c

          lhs[,malo] <- xc
          medidashs[malo] <- fxc


        }else{                 # Shrink

          # Paso 6: Shrink 
          #######################################################
          # we shrink the simplex fixing the best candidate.
          # Proportion: sigma.

          xbueno <- lhs[,bueno]
          lhs[,-bueno] <- apply(lhs[,-bueno], 2, function(x) (x + xbueno)/2)

          # for(i in (1:ncol(lhs))[-bueno]){
          #   cat("Probando h=",exp(lhs[,i])+htol,"\n") ###
          # }

          medidashs[-bueno] <- apply(exp(lhs[,-bueno])+htol, 2, function(h) medida(h))

          #- end: shrink -#######################################

        }


        #- end: reflection -#########################################################



        
      }

      # ----------------------------------------------------------------------------

    }else{ # f(h_r) < f(h_feo)


      # if h_r is better than the best, h_bueno, replace the later by the former.

      if(fxr >= fbueno){

        lhs[,malo] <- xr
        medidashs[malo] <- fxr

      }else{   # in case f(h_r) > f(h_bueno) ==> Step 4: expansion                     


        # Step 4: expansion
        ##################################################3
        # The search is extended further in the direction of
        #  h_r from h_0, computing a new candidate h_e,
        #
        #            h_e <················ h_r ----------------- h_0
        # Proportion: |------------------ gamma ------------------|

        xe <- 2*xr - x0              # Expansion

        # cat("Probando h=",exp(xe)+htol,"\n") ###

        fxe <- medida(exp(xe)+htol) # the measurement is computed.


        if(fxe < fxr){ # if the expanded candidate, h_e, is better than the reflection, h_r
                       # then replace the worst candidate, h_malo, with h_e.

          lhs[,malo] <- xe
          medidashs[malo] <- fxe

        }else{ # if the expanded candidate, h_e, is worst than the reflection, h_r
               # then replace the worst candidate, h_malo, with h_r. 
          lhs[,malo] <- xr
          medidashs[malo] <- fxr
        }

      }
      


    }

    #- end: check ---------------------------------------------------------------------

    # plot(t(exp(lhs)+htol), pch = 19, col = "red", xlim=c(0,h0[1]+3*as),ylim=c(0,h0[2]+3*as)) ###
    # lines(t(exp(cbind(lhs,lhs[,1]))+htol), col = "red") ###


    it=it+1 # the iteration counter is increased.

    lhb <- lhs[,which.min(medidashs)]

    if(it == maxiter) break # if the maximum allowed number of iteration is reached, end.


    if(sqrt(sum(((exp(lhb)-exp(lhb0))/exp(lhb))^2)) < eps2) j <- j+1 else j <- 0  # si mejoramos, reiniciamos el contador

    if(j == maxiter.tol) break # if there is not improvement in the pre-set tolerance number of iteration, ends.
    
    lhb0 <- lhb

    # cat("Desv. Típica=",sd(medidashs),", it=",it,", j=",j,"\n") ###

  }
  return(exp(lhb0)+htol)
  # return(exp(lhs[,which.min(medidashs)])+htol)
}

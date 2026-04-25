PMS1 <- function(X,Y,x,ymalla,h1,h2,eps,k,n){
 .Call("PMF1",as.numeric(X),as.numeric(Y)
                               ,as.numeric(x), as.numeric(ymalla),as.numeric(h1),as.numeric(h2)
                               ,as.numeric(eps),as.integer(k),as.integer(n))
}

norma <- function(v,V){
  return(sqrt(apply((v-V)^2,1,sum)))
}
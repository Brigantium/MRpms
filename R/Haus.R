.Haus<-function(A,B){max(apply(abs((outer(A,B, "-"))),2,min))}

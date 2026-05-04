set.seed(2026)
x1 <- 10*runif(100)
m1 <- sin(x1) + 0.3*rnorm(100)
x2 <- 10*runif(100)
m2 <- sin(x2) + 0.3*rnorm(100) - 3
X <- c(x1, x2)
Y <- c(m1, m2)
Ejemplo1 <- cbind(X, Y)

usethis::use_data(Ejemplo1, overwrite = TRUE)

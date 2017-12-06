#' @title Right-Truncated Poisson Probability Mass Function
#' 
#' @description PMF of a right truncated Poisson distribution.
#' 
#' @param y The vector of counts.
#' @param lambda The vector of means
#' @param c The cutoff of inflation
#' 
#' @author Michael Floren
dtpois <- function(y, lambda, c){
  n <- length(y)
  out <- numeric(length=n)
  b <- numeric()
  for(i in 1:n){
    for(j in 0:c){
      b[j+1] <- ((exp(-lambda[i])*lambda[i]^j)/(factorial(j)))
    }
    b <- sum(b)
    out[i] = dpois(x=y[i], lambda=lambda[i]) / b
  }
  out[y>c] <- 0 #this PMF is 0 by definition for y>c
  out
}
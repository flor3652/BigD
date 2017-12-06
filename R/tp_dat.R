#' @title Generation of truncated Poisson data
#' 
#' @description This function generates truncated Poisson data, with a truncation at c (counts can include c).
#' 
#' @param n The number of data points to be generated.
#' @param lambda The mean parameter for the truncated Poisson.
#' @param c The cutoff of inflation.
#' 
#' @details This generates data from a right truncated Poisson distribution, where generated y's are in the range \eqn{0 \le y \le c}.
#' 
#' @examples 
#' #show a histogram of right truncated Poisson data
#' hist(tp_dat(n=3000, lambda=3, c=5))
#' 
#' @author Michael Floren
#' @md

tp_dat <- function(n, lambda, c){ #truncated Poisson
  if(n != length(lambda) & length(lambda)>1) 
    stop("The length of lambda must be 1 or match n.") #only support if a different lambda is given for each n, or the same lambda is given for all of them...
  if(c <= lambda)
    stop("You're being an idiot (lambda is greater than c)") #doesn't make sense for the mean of the distribution to be larger than the cutoff (tell user their being an idiot)
  
  if(length(lambda)!=n){ #this would be the length(lambda)==1 and n>1 case...
    lambda <- rep(lambda, n)
  }
  
  out <- numeric()
  for(j in 1:length(lambda)){
    probs <- numeric()
    for(i in 0:c){
      probs[i+1] <- ((exp(-lambda[j])*lambda[j]^i)/(factorial(i)))
    }
    probs = probs/(sum(probs))
    #sample(seq(0,c), n, replace=TRUE, prob=probs)
    out[j] <- sample(0:c, size=1, prob=probs)
  }
  
  out
}

#' @title Generation of Small Count Inflated Poisson Data
#' 
#' @description This function generates small count inflated Poisson (SCIP) data, giving a list of values and group membership.
#' 
#' @param n The number of data points to be generated.
#' @param desmat The design matrix to use for data generation.
#' @param sam A vector of parameters to use for the small count distribution.
#' @param kop A vector of parameters to use for the large count distribution.
#' @param gam A vector of parameters to use for the binomial distribution (1 is small count group).
#' @param c The cutoff of inflation
#' 
#' @details This function takes a set of arguments and returns a vector of two elements. The first element is the generated count from the SCIP distribution, the second element is the group of membership (for dissertation purposes). The returning of the group membership element may be removed at a later date.
#' 
#' @return This function currently returns a list of two vectors
#'   \item{y}{A vector of the generated count variables}
#'   \item{g}{A vector of the group each generated count belongs to}
#'   
#' @examples
#' #show a histogram of small count inflated data
#' hist(scip_dat(n=3000, desmat=matrix(rep(1, 3000), ncol=1), sampi=log(.1), koppa=log(4), gamma=log(1), c=3)$y)
#' 
#' @author Michael Floren

scip_dat <- function(n, desmat, sam, kop, gam, c){
  if(!all.equal(ncol(desmat), length(sam), length(kop), length(gam)))
    stop("Dimensional issues between design matrix, sampi, koppa, and gamma. Please check dimensions (ncol of design matrix should match length of sampi, koppa, and gamma).")
  
  lam_s <- exp(desmat%*%sam) #lambda_s
  lam_l <- exp(desmat%*%kop) #lambda_l
  pi <- exp(desmat%*%gam) / (1+exp(desmat%*%gam)) #pi_s
  
  obs <- numeric()
  group <- numeric()
  for(i in 1:n){
    group[i] <- rbinom(1,1,pi)
    if(group[i] == 1){
      obs[i] <- tp_dat(n=1, lambda=lam_s[i], c=c)
    } else{
      obs[i] <- rpois(n=1, lambda=lam_l[i])
    }
  }
  list(y=obs, g=group)
}

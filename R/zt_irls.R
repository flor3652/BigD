#' @title The IRLS Function
#' 
#' @description
#' This function is designed to perform an IRLS algorithm to calculate parameters for regression on a given dataset (this will be designed to handle count data).
#' 
#' @param formula A formula of the equation that should be used for regression.
#' @param data The data that should be used. Formula should use variable names as seen in the data.
#' @param epsilon This is the criterion for convergence: set manually, for getting a feel for how it works
#' @param maxit The maximum number of iterations you are willing to sit through before you blow your brains out...
#' @param ... Not currently supported.
#' 
#' @details
#' This is a practice function to see if I'm understanding the concepts (and if it actually runs). Documentation is provided purely for practice.
#' 
#  For whatever reason, I CANNOT get this link to the "rm" function to actually create a hyperlink
#' @seealso \code{\link[pscl]{zeroinfl}}
#' 
#' @examples 
#' # using the "fish" dataset, downloaded with this package
#' # ztfish <- read.csv("http://www.ats.ucla.edu/stat/data/fish.csv") #if you want to grab it manually
#' # ztfish <- ztfish[ztfish$count>0,]
#' 
#' #an example with a two way interaction
#' irls(count~persons*camper, ztfish)
#' 
#' #an example with a three way interaction
#' irls(count~persons*camper*nofish, ztfish)
#' 
#' @author Michael Floren

zt_irls <- function(formula, data, epsilon = .001, maxit=10, ...){
  ### Step 1 (steps are from the paper)
  initmod <- glm(formula, data=data, family="poisson")
  initcoef <- initmod$coefficients
  alpha <- matrix(nrow=length(initcoef)) # this is the coefficient matrix as listed in the reference article (see the practice script). Each iteration of coefficiencts is a new column: to match what the article calls  "alpha" (for the matrix multiplication equations, we don't have to artifically transpose for formatting)...
  rownames(alpha) <- names(initcoef)
  alpha[,1] <- initcoef #since we already know these, go on and enter them as the first iteration
  
  ### Step 2 (I think the second half of this should be inside of a loop (after z is fully defined))
  z <- cbind(intercept=1, as.matrix(initmod$model[,-1]))
  if(sum(attributes(initmod$terms)$order>1)){
    intterms <- colnames(attributes(initmod$terms)$factors)[attributes(initmod$terms)$order>1]
    int <- matrix(ncol=length(intterms), nrow=nrow(z))
    for(i in seq(1,length(intterms))){
      int[,i] <- apply(X = z[,colnames(z) %in% stringr::str_split(intterms,":")[[i]]], MARGIN = 1, prod)
    }
    colnames(int) <- intterms
    z <- cbind(z,int)
  }
  
  k <- 1 #iterator to compare to maxit
  l <- numeric()
  l[1] <- epsilon + 1 #arbitrary start...
  while(k<maxit){
    loz <- exp(z %*% alpha[,k, drop=FALSE]) #lambda of z... for the iterations, the i^th column.
    mu <- loz/(1-exp(-loz))
    D <- matrix(0,ncol=length(loz), nrow=length(loz))
    diag(D) <- (loz*(1-exp(-loz)*(1+loz)))/((1-exp(-loz))^2)
    Dinv <- solve(D)
    t <- z %*% alpha[,k, drop=FALSE] + Dinv %*% (initmod$y-mu)#Assume that Z is z (I guess?). This doesn't seem right, as the previous z in the paper was transposed befor multiplied with alpha, but thats what I'll do!
    l[k+1] <- t(t - z %*% alpha[,k, drop=FALSE]) %*% D %*% (t - z %*% alpha[,k, drop=FALSE])
    
    # Step 3 ----
    alpha = cbind(alpha,solve(t(z)%*%D%*%z) %*% t(z) %*% D %*% t) #adding on another column (of the updated coefficients)...
    k <- k+1
    if(abs(l[k]-l[k-1])/l[k-1]<epsilon){
      print(paste("Convergence within epsilon of", epsilon, "reached after", k, "iterations."))
      break
    }
  }
  if(k>maxit)
    print(paste("Convergence not reached after", k, "iterations (the maximum set by user)."))
  
  out <- list(coefficients=alpha[,k], iterations=k) #could also do alpha to see all of the iteracted coefficients...
  out
}



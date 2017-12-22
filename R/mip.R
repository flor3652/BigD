#' @title Multinomially Inflated Poisson Model
#' 
#' @description This is the multinomially inflated Poisson distribution as demonstrated by Giles (2007).
#' 

mip <- function(y, desmat, c){
  # log-liklihood given outcome, parameters, and design matrix...
  ll <- function(y, desmat, param){
    #reading in gammas and betas from the long list of parameters
    p <- ncol(desmat)
    c <- length(param)/ncol(desmat) - 2
    param_mat <- matrix(ncol=c+3, nrow=p)
    for(i in 1:(c+1)){ #so this is 0:c (doesn't do gamma_J)
      param_mat[,i] <- param[(i*p-(p-1)):(i*p)]
    }
    param_mat[,c+2] <- rep(0,p)#this is gamma_J
    param_mat[,c+3] <- param[((c+2)*p-(p-1)):((c+2)*p)] 
    colnames(param_mat) <- c(paste0("gamma",0:(c+1)), "beta") #remember: gamma_(c+1) is gamma_J (0's)
    
    # prepping omega calculation
    denom_for_omega_mat <- matrix(ncol=c+1, nrow=length(y))
    for(i in 0:c){ #for all the gammas except for gamma zero, cause this is what the denom sum is. gamma_(c+1) is gamma_J (so including J). Gamma_J is zero.
      denom_for_omega_mat[,i+1] <- exp(desmat %*% param_mat[,paste0("gamma",i)])
    }
    denom_for_omega <- 1 + apply(denom_for_omega_mat, 1, sum)
    
    # calculating omegas : something weird going on here (matches the long hand through the log-likelihood output) ******************************************
    omega_mat <- matrix(ncol=c+2, nrow=length(y))
    for(i in 0:(c+1)){ #creating an omega_(c+1) (aka, omega_J)
      omega_mat[,i+1] <- exp(desmat %*% param_mat[,paste0("gamma",i)])/denom_for_omega #right now, desmat %*% param (not exponentiated) looks to be a good prediction of the weights (for whatever reason) (if you use the optimized values for your parameters and run through this line by line...). Almost just proportional: cbind(desmat[,2],desmat%*%param_mat[,1], desmat[,2]-desmat%*%param_mat[,1], (desmat[,2]-desmat%*%param_mat[,1])/desmat%*%param_mat[,1])
    }
    colnames(omega_mat) <- paste0("omega", 0:(c+1)) #remember: omega_(c+1) = omega_J
    
    #calculating coefficient for poisson 
    coef_p <- 1-apply(omega_mat[,-ncol(omega_mat), drop=FALSE], 1, sum)# take off the last column (the c+2nd column (omega_(c+1) aka omega_J)) from omega mat
    
    #calculating pois
    lambda <- exp(desmat %*% param_mat[,"beta"])
    pois <- (exp(-lambda)*lambda^y)/(factorial(y))
    
    ll_out_vec <- numeric(length=c+2) #this is each of the sums in the log likelihood from Giles (0 to J (which is c+1)), listed as a vector (eventually will be summed for the log-likelihood).
    for(i in 0:c){
      ll_out_vec[i+1] <- sum((y==i)*log((omega_mat[,paste0("omega",i)] + (coef_p)*pois))) #add everything up (add 0's unless y is i)
    }
    ll_out_vec[c+2] <- sum((y>c)*log((coef_p)*pois)) #doing the poisson alone by hand
    ll_out <- sum(ll_out_vec) #matches with longhand...
    ll_out
  }
  
  p <- ncol(desmat)
  
  #starting guess for parameters (just use 0 vector I guess?)
  initial_guess_param <- rep(0, p*(c+2)) #gammas from 0 through c, then beta (each has p parameters)...
  
  of <- function(params){ #optimization function (just a function of parameters)
    ll(y=y, desmat=desmat, param=params)
  }
  
  full_out <- optim(initial_guess_param, of, control = list(fnscale=-1))
  est <- full_out$par
  convergence_issue <- full_out$convergence!=0
  
  
  out <- list()
  for(i in 1:(c+1)){ #so this is gamma_0:c
    out[[i]] <- est[(i*p-(p-1)):(i*p)]
    names(out)[i] <- paste0("logistic",i-1)
  }
  out[[c+2]] <- est[((c+2)*p-(p-1)):((c+2)*p)]
  names(out)[c+2] <- "beta"
  out$convergence_issue <- convergence_issue
  out
}
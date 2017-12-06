#' @title SCIP Function
#' 
#' @description Fits a Small Count Inflated Poisson model.
#' 
#' @author Michael Floren

scip <- function(y, desmat, c, conv=1e-9, maxit=100, progressBar=TRUE, timeUnit="auto"){
  ### Functions used later on ###
  tau_s <- function(y, desmat, param, c){ #take the design matrix with sampi, koppa, gamma, and the design matrix (written short as "gamma" is a system word)
    sam <- param[1:ncol(desmat)]
    kop <- param[(ncol(desmat)+1):(2*ncol(desmat))]
    gam <- param[(2*ncol(desmat)+1):(3*ncol(desmat))]
    lam_s <- exp(desmat %*% sam)
    lam_l <- exp(desmat %*% kop)
    pi <- exp(desmat%*%gam) / (1+exp(desmat%*%gam))
    scip_pmf <- pi*dtpois(y=y, lambda=lam_s, c=c) + (1-pi)*dpois(x=y, lambda=lam_l)
    
    out <- (pi * dtpois(y=y, lambda=lam_s, c=c)) / scip_pmf
    out[y>c] <- 0 #if the outcome is greater than c, manually set the tau to be zero (it should be anyways, but the auto-zero of the dtpois has been removed for problems in the logarithm)
    out
  }
  
  Q <- function(y, desmat, tau, param, c){
    sam <- param[1:ncol(desmat)]
    kop <- param[(ncol(desmat)+1):(2*ncol(desmat))]
    gam <- param[(2*ncol(desmat)+1):(3*ncol(desmat))]
    lam_s <- exp(desmat %*% sam)
    lam_l <- exp(desmat %*% kop)
    pi <- exp(desmat%*%gam) / (1+exp(desmat%*%gam))
    
    sum(log((pi*dtpois(y=y, lambda = lam_s, c = c))^tau) + 
          log(((1-pi)*dpois(x=y, lambda = lam_l))^(1-tau)))
  }
  
  log_pdf_scip <- function(y, desmat, tau, param, c){
    sam <- param[1:ncol(desmat)]
    kop <- param[(ncol(desmat)+1):(2*ncol(desmat))]
    gam <- param[(2*ncol(desmat)+1):(3*ncol(desmat))]
    lam_s <- exp(desmat %*% sam)
    lam_l <- exp(desmat %*% kop)
    pi <- exp(desmat%*%gam) / (1+exp(desmat%*%gam))
    pois <- dpois(x=y, lambda=lam_l) #so that you only calculate it once
    
    sum(log(pi*(dtpois(y=y, lambda=lam_s, c=c) - pois) + pois))
  }
  
  
  
  ### Actual implementation ###
  # Creating estimates for sampi, koppa, and gamma. Tracking all iterations currently (don't have to do it this way: could just wipe the previous iteration out on each run...)
  est_param <- matrix(ncol=3*ncol(desmat))
  init_est_sampi <- coef(glm(y~desmat[,-1], family=poisson))
  init_est_koppa <- coef(glm(y~desmat[,-1], family=poisson))
  init_est_gamma <- coef(glm(as.numeric(y<c)~desmat[,-1], family=binomial))
  est_param[1,] <- c(init_est_sampi, init_est_koppa, init_est_gamma)
  
  # Creating initial weights: not tracking these for each iteration
  tau <- tau_s(y=y, desmat=desmat, param=est_param[1,], c=c)
  
  # Setting the start time (eventually can give a number of how long this took)
  starttime <- Sys.time()
  
  # Use progress bar version if a progress bar was requested: pb is out of the max number of iterations
  if(progressBar) {
    cat(paste0("Percentage of maximum iterations completed (i=",maxit,"):\n"))
    pb <- txtProgressBar(style=3)
  }
    
  # Run EM
  timing <- numeric() #for large data. take this and tic-toc and timing message out for regular runs...
  for(i in 2:maxit){
    tic <- Sys.time()
    #optimization function for parameters (only takes parameters as arguments)
    of_param <- function(param){
      Q(y=y, desmat=desmat, tau=tau, param=param, c=c)
    }
    #solving the derivative for the current iteration
    est_param <- rbind(est_param, optim(est_param[i-1,], of_param, control = list(fnscale=-1))$par)
    
    #updating tau
    tau <- tau_s(y=y, desmat=desmat, param = est_param[i,], c=c)
    
    #checking if convergence has been reached
    check <- max(abs(est_param[i,] - est_param[i-1,]))
    if(check < conv){
      break
    }
    if(progressBar)
      setTxtProgressBar(pb,value=i/maxit)
    toc <- Sys.time()
    timing[i-1] <- difftime(toc,tic, units="s")
    cat("Finished",i,"th iteration at", Sys.Date(),". This iteration took", timing[i-1],"seconds. Average run-time per iteration is", mean(timing), "seconds.\n")
  }
  if(progressBar){
    close(pb)
    if(i == maxit){
      message(paste0("Convergence not met after maximum number of iterations (i=",maxit,")"))
    } else {
      message("Convergence met after ", i, " iterations.")
    }
  }
  endtime <- Sys.time()
  iterations <- i
  
  # setting output
  est_sam <- est_param[i, 1:ncol(desmat)]
  est_kop <- est_param[i, (ncol(desmat)+1):(2*ncol(desmat))]
  est_gam <- est_param[i, (2*ncol(desmat)+1):(3*ncol(desmat))]
  runtime <- difftime(endtime, starttime, units = timeUnit) #add |, units='s'| to force this into seconds...
  
  # out
  list(param=est_param, sc = est_sam, lc=est_kop, logistic=est_gam, runtime=runtime, iterations=iterations)
}



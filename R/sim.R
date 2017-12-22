#' @title Simulation Function for Dissertation
#' 
#' @description This function will run the simulation for my dissertation. Some of the aspects of the function are specifically designed for this dissertation (such as population parameters), and won't be included as arguments.
#' 
#' 
#' 

sim <- function(n, c, k=1000, meta_time_unit="s", progress=TRUE){
  p <- 2
  mean_of_IV <- 197.7828
  sd_of_IV <- 40.21599
  
  # setting different true values depending on c
  if(c==2){
    true_sampi <- c(-1.4890070090395, 0.000896026713969103)
    true_koppa <- c(2.26780053964254, -0.000373799698134193)
    true_gamma <- c(1.79118900609313, -0.000338077141916478)
  } else if (c==6){
    true_sampi <- c(-0.991816836348344, 0.000632613090524932)
    true_koppa <- c(2.51780279408921, -0.000375747468669122)
    true_gamma <- c(2.1998559804243, -0.00025944850892789)
  } else if (c==8){
    true_sampi <- c(-0.989457553026938, 0.000623795679864685)
    true_koppa <- c(2.51994084026393, -0.000385423153504593)
    true_gamma <- c(2.20836925037315, -0.000300453184670207)
  } else{
    stop(paste0("Population parameters for c=",c," not set. Simulation stopped."))
  }
  
  
  
  #creating the place to save all of the results. These will have elements named for each piece of the output (e.g., sc, lc, logistic), where each row is a different run
  zip_results <- list()
  scip_results <- list()
  mip_results <- list()
  meta <- list()
  
  #creating list elements for each piece of output
  for(i in 1:2){
    zip_results[[i]] <- matrix(ncol=p, nrow=k) #columns is the number of parameters (in this case, intercept and slope)
  }
  for(i in 1:3){
    scip_results[[i]] <- matrix(ncol=p, nrow=k)
  }
  for(i in 1:(c+2)){ #0-c and beta
    mip_results[[i]] <- matrix(ncol=p, nrow=k)
  }
  names(zip_results) <- c("count", "logistic") #the names and order from zeroinfl
  names(scip_results) <- c("sc","lc","logistic")
  names(mip_results) <- c(paste0("logistic",0:c),"beta")
  
  #creating a place to store meta-information (iterations for SCIP, runtime for all) for each run and overall. NA's by default (don't want 0's by default, as they may not get noticed as errors (if they aren't overwritten, for some reason))
  zip_results$runtime <- as.numeric(rep(NA,k))
  zip_results$convergence_issue <- as.numeric(rep(NA,k))
  mip_results$runtime <- as.numeric(rep(NA,k))
  mip_results$convergence_issue <- as.numeric(rep(NA,k))
  scip_results$runtime <- as.numeric(rep(NA,k))
  scip_results$iterations <- as.numeric(rep(NA,k))
  scip_results$convergence_issue <- as.numeric(rep(NA,k))
  meta$runtime <- as.numeric(rep(NA,k))
  
  #creating an error log for each of the methods
  zip_results$error_log <- data.frame(error_bin=rep(0,k), error_msg=character(length=k), stringsAsFactors=FALSE)
  mip_results$error_log <- data.frame(error_bin=rep(0,k), error_msg=character(length=k), stringsAsFactors=FALSE)
  scip_results$error_log <- data.frame(error_bin=rep(0,k), error_msg=character(length=k), stringsAsFactors=FALSE)
  
  #performing the trials
  for(i in 1:k){
    #starting the run clock
    meta_starttime <- Sys.time()
    
    
    
    # generating the dsign matrix
    desmat <- cbind(rep(1,n), rnorm(n, mean=mean_of_IV, sd=sd_of_IV)) #these are the population parameters from weight in the BRFSS (agam, manually entered)
    
    # generating data
    y <- scip_dat(n=n, desmat=desmat, sam=true_sampi, kop=true_koppa, gam=true_gamma, c=c)$y
    
    
    
    ### Fitting the ZIP ####
    # Making a dataset for the ZIP 
    zip_dat <- cbind(y=y,x=desmat[,-1])
    
    # fitting the zip model
    #library(pscl) #eventually need to take this out
    zip_start<-Sys.time()
    zip_attempt <- tryCatch(zip_mod <- pscl::zeroinfl(y~., data=as.data.frame(zip_dat)), error=function(e) c(error_bin=1,e))
    zip_stop<-Sys.time()
    
    if("error_bin" %in% names(zip_attempt)){ #if an error
      zip_results$error_log$error_bin[i] <- zip_attempt$error_bin
      zip_results$error_log$error_msg[i] <- zip_attempt$message
    } else { #record the fit
      zip_results$count[i,] <- zip_mod$coefficients$count
      zip_results$logistic[i,] <- zip_mod$coefficients$zero
      zip_results$runtime[i] <- difftime(zip_stop, zip_start, units=meta_time_unit)
      zip_results$convergence_issue[i] <- 1-zip_mod$converged
    }
    
    
    
    ### Fitting the SCIP ####
    scip_start<-Sys.time()
    scip_attempt <- tryCatch(scip_mod <- scip(y=y, desmat=desmat, c=c), error=function(e) c(error_bin=1,e))
    scip_stop<-Sys.time()
    
    #error checking and recording results (if no error)
    if("error_bin" %in% names(scip_attempt)){ #if an error
      scip_results$error_log$error_bin[i] <- scip_attempt$error_bin
      scip_results$error_log$error_msg[i] <- scip_attempt$message
    } else { #record the fit
      scip_results$sc[i,] <- scip_mod$sc
      scip_results$lc[i,] <- scip_mod$lc
      scip_results$logistic[i,] <- scip_mod$logistic
      scip_results$runtime[i] <- difftime(scip_stop, scip_start, units=meta_time_unit)
      scip_results$iterations[i] <- scip_mod$iterations
      scip_results$convergence_issue[i] <- scip_mod$convergence_issue
    }
    
    
    ### Fitting the MIP ####
    mip_start<-Sys.time()
    mip_attempt <- tryCatch(mip_mod <- mip(y=y, desmat=desmat, c=c), error=function(e) c(error_bin=1,e))
    mip_stop<-Sys.time()
    
    #error checking and recording results (if no error)
    if("error_bin" %in% names(mip_attempt)){ #if an error
      mip_results$error_log$error_bin[i] <- mip_attempt$error_bin
      mip_results$error_log$error_msg[i] <- mip_attempt$message
    } else { #record the fit
      # as they should all be named the same, using the names from one to iterate over the other...
      for(j in c(paste0("logistic",0:c),"beta")){ #the names of the mip_results columns that don't have to do with meta stuff...
        mip_results[[j]][i,] <- mip_mod[[j]]
      }
      mip_results$runtime[i] <- difftime(mip_stop, mip_start, units=meta_time_unit)
      mip_results$convergence_issue[i] <- mip_mod$convergence_issue
    }
    
    
    
    #ending the run clock
    meta_endtime <- Sys.time()
    
    #Recording meta-time (time for the full run)
    meta$runtime[i] <- difftime(meta_endtime, meta_starttime, units = meta_time_unit)
    
    if(progress)
      cat(paste0("\r",round(i/k*100),"% complete (working on ", toOrdinal::toOrdinal(i+1)," run). Last run took ", round(meta$runtime[i],2), " ", meta_time_unit,". Average run is ", round(mean(meta$runtime, na.rm=TRUE),2), " ", meta_time_unit, ". Predicted completion time is ", Sys.time()+round(mean(meta$runtime, na.rm=TRUE))*(k-i),"."))
  }
  
  
  
  list(zip_results=zip_results, scip_results=scip_results, mip_results=mip_results, meta=meta, n=n, c=c)
}

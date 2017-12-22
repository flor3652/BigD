#' @title Prediction Function for Dissertation
#' 
#' 
#' @description This function will run the prediction for my dissertation. The goal is to establish prediction accuracy for each of the models. So, we'll need to generat a design matrix, generate some data, have a predicted value based on the design matrix, have an actual value, then also have an actual group vs a predicted group for each of the models...
#' 
#' I want to matrix multiple all of the things, then use logistic to determine which one I should actually use...
#' 
#' @param n The sample size to use. This should be 20\% of the actual data used in the simulation.
#' @param c The cutoff to use. This should match that used in the simulation.
#' @param sim_dat Data from the simulation function.
#' 
#' @author Michael Floren
pred <- function(sim_dat){
  ### initial parameters ###
  percent_for_test <- .2 #test sample size is 20% of the actual sample. May choose to use 25% instead, but not a big deal...
  n <- sim_dat$n*percent_for_test
  c <- sim_dat$c
  
  ### functions ###
  mse <- function(y,yhat){ #mean squared error
    sum((yhat-y)^2)/length(y)
  }
  
  #the next couple functions need a threshhold to compare to
  thresh <- .5 
  
  pc <- function(g, ghat){ #percent correct, given decimals for ghat. Remember: 1 is small count, 0 is large count (this works either way, but g and ghat have to be consistent...)
    sum((ghat>thresh)==g)/length(g)
  }
  
  auc <- function(g, ghat){
    pROC::auc(pROC::roc(as.numeric(g),as.numeric(ghat)))
  }
  
  # all overall numbers should be the same as for the simulation
  k <- nrow(sim_dat$zip_results$count) #the number of rows for any of these should be the same...
  p <- ncol(sim_dat$zip_results$count)
  mean_of_IV <- 197.7828 
  sd_of_IV <- 40.21599
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
  
  # making generic lists to hold the information
  verbose_results <- list()
  results <- list()
  
  #defining matrices for the results
  for(i in 1:3){
    results[[i]] <- as.data.frame(matrix(ncol=3, nrow=k))
    colnames(results[[i]]) <- c("zip","mip","scip")
  }
  names(results) <- c("mse","percent_cor", "auc")
  results$weighted_mse <- as.data.frame(matrix(ncol=3, nrow=k))
  colnames(results$weighted_mse) <- c("zip", "scip", "mip")
  
  for(i in 1:k){
    # generating data for the run
    desmat <- cbind(rep(1,n), rnorm(n, mean=mean_of_IV, sd=sd_of_IV))
    actual <- scip_dat(n=n, desmat=desmat, sam=true_sampi, kop=true_koppa, gam=true_gamma, c=c)
    
    verbose_results[[i]] <- as.data.frame(matrix(ncol=p+8, nrow=n)) #the verbose results for the iteration/run. Columns: p for desmat, 2 for out (group and count), 2 for predictions (count and group) for ZIP, SCIP, and MIP
    colnames(verbose_results[[i]]) <- c(paste0("desmat", seq(0,ncol(desmat)-1)), "y", "g", "yhat_zip", "ghat_zip", "yhat_scip", "ghat_scip", "yhat_mip", "ghat_mip")
    
    #setting info for the run
    for(j in 1:ncol(desmat))
      verbose_results[[i]][,paste0("desmat",j-1)] <- desmat[,j]
    
    verbose_results[[i]]$y <- actual$y
    verbose_results[[i]]$g <- actual$g
    
    
    
    #### Predictions: ZIP ####
    zip_group_pred <- exp(desmat%*%sim_dat$zip_results$logistic[i,])/(1+exp(desmat%*%sim_dat$zip_results$logistic[i,])) #exp(x\beta)/(1+exp(x\beta))
    zip_count_pred <- exp(desmat%*%sim_dat$zip_results$count[i,])
    
    #setting the weighted mse
    results$weighted_mse$zip[i] <- mse(y=actual$y, 
                                       yhat=zip_group_pred*0 + (1-zip_group_pred)*zip_count_pred)
    
    #for logistic results that point towards zero, set the predicted count to 0 (use same threshold as above)
    zip_count_pred[zip_group_pred>thresh] <- 0
    
    #set the verbose results information
    verbose_results[[i]]$yhat_zip <- zip_count_pred
    verbose_results[[i]]$ghat_zip <- zip_group_pred
    
    #set the results information
    results$mse$zip[i] <- mse(y=actual$y, yhat=zip_count_pred)
    results$percent_cor$zip[i] <- pc(g=actual$g, ghat=zip_group_pred)
    results$auc$zip[i] <- auc(g=actual$g, ghat=zip_group_pred)
    
    
    
    #### Predictions: SCIP ####
    # Just treating the logistic piece as the group prediction...
    scip_group_pred <- exp(desmat%*%sim_dat$scip_results$logistic[i,])/(1+exp(desmat%*%sim_dat$scip_results$logistic[i,])) #the logistic piece
    scip_sc_pred <- exp(desmat%*%sim_dat$scip_results$sc[i,])
    scip_lc_pred <- exp(desmat%*%sim_dat$scip_results$lc[i,])
    scip_count_pred <- ifelse(scip_group_pred>thresh, scip_sc_pred, scip_lc_pred)
    
    #set the verbose results information
    verbose_results[[i]]$yhat_scip <- scip_count_pred
    verbose_results[[i]]$ghat_scip <- scip_group_pred
    
    #set the results information
    results$mse$scip[i] <- mse(y=actual$y, yhat=scip_count_pred)
    results$percent_cor$scip[i] <- pc(g=actual$g, ghat=scip_group_pred)
    results$auc$scip[i] <- auc(g=actual$g, ghat=scip_group_pred)
    results$weighted_mse$scip[i] <- mse(y=actual$y,
                                        yhat= scip_group_pred*scip_sc_pred + (1-scip_group_pred)*scip_lc_pred)
    
    
    
    #### Predictions: MIP ####
    names_of_logistic_pieces <- paste0("logistic", 0:c)
    
    #doing the denominator first
    denom_for_logistic_mat <- matrix(ncol=c+1, nrow=length(y))
    for(j in 1:(c+1)){ #for all the gammas except for gamma zero, cause this is what the denom sum is. gamma_(c+1) is gamma_J (so including J). Gamma_J is zero.
      denom_for_logistic_mat[,j] <- exp(desmat %*% sim_dat$mip_results[[names_of_logistic_pieces[j]]][i,])
    }
    denom_for_logistic <- 1 + apply(denom_for_logistic_mat, 1, sum)
    
    #determining the logistic prediction pieces
    mip_logistic <- as.data.frame(matrix(ncol=c+2, nrow=n))
    colnames(mip_logistic) <- c(names_of_logistic_pieces, "logisticJ")
    for(j in 1:(c+1)){
      mip_logistic[,j] <- exp(desmat%*%sim_dat$mip_results[[names_of_logistic_pieces[j]]][i,])/denom_for_logistic
    }
    mip_logistic$logisticJ <- 1-apply(mip_logistic[,-ncol(mip_logistic)], 1, sum)
    mip_group_pred <- 1-mip_logistic$logisticJ #the probability of being in the small count is 1-\pi_l
    mip_count_pred <- exp(desmat%*%sim_dat$mip_results$beta[i,])
    
    #setting the weighted MSE information
    weighted_mip <- mip_logistic
    for(j in 1:(c+1))
      weighted_mip[,j] <- mip_logistic[,j]*(j-1) #the column times its count
    weighted_mip$logisticJ <- mip_logistic$logisticJ * mip_count_pred #the count probability times the count prediction
    results$weighted_mse$mip[i] <- mse(y=actual$y,
                                       yhat=apply(weighted_mip, 1, sum))
    
    #for each row, check which logistic piece is higher
    group <- apply(mip_logistic, 1, function(x) which(x==max(x))[1])-1 #if multiple are tied, just take the first one... subtract 1 to match the column index with the meaning (first column is a count of zero)
    for(j in 0:c){ #don't do J, as if the max group is J we want to leave it alone. Only change for 
      mip_count_pred[group==j] <- j
    }
    
    #set the verbose results information
    verbose_results[[i]]$yhat_mip <- mip_count_pred
    verbose_results[[i]]$ghat_mip <- mip_group_pred
    
    #set the results information
    results$mse$mip[i] <- mse(y=actual$y, yhat=mip_count_pred)
    results$percent_cor$mip[i] <- pc(g=actual$g, ghat=mip_group_pred)
    results$auc$mip[i] <- auc(g=actual$g, ghat=mip_group_pred)
  }
  
  list(verbose_results=verbose_results, results=results)
}
# Common Beta version of the Weibull MCMC

library(coda)
library(Rcpp)

cppFunction('double logWeibullSumCPP(NumericVector x, double shape, double scale) {
            
            int nx = x.size();
            double d = 0.0;
            
            for(int i = 0; i < nx; i++) {
            d = d + log(shape) + log(scale) + (shape - 1)*log(x[i]) - scale*pow(x[i],shape);
            }
            
            return d;
            }')


################################################################################
# Start of MCMC function
################################################################################

WeibullOneShapeMCMC <- function(data, samples = 500, shapePriorA = .001,
                                shapePriorB = .001, HyperG1 = .001, HyperG2 = .001,
                                hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                                hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                                alphaStart = 1, betaStart = 1, theta1Start = 1,
                                theta2Start = 1, shapeStart = 1, tuningA = 1,
                                tuningS = 1, trun = T, EX = F){
  
  tuningAV <- rep(tuningA, 26)
  
  # log posterior functions for MCMC draws
  logPostComp <- function(d, parm, shape) {
    
    lp <- 0
    
    for(k in 5:length(parm)){
      lp <- lp + 
        logWeibullSumCPP(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 0 & d[[k-4]]$Phase == 1 & d[[k-4]]$trun == F)], shape, parm[k]) +
        logWeibullSumCPP(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 0 & d[[k-4]]$Phase == 2 & d[[k-4]]$trun == F)], shape, parm[k]*parm[3]) +
        logWeibullSumCPP(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 0 & d[[k-4]]$Phase == 3 & d[[k-4]]$trun == F)], shape, parm[k]*parm[3]*parm[4]) -
        sum(parm[k]*(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 1 & d[[k-4]]$Phase == 1 & d[[k-4]]$trun == F)]^shape)) -
        sum(parm[k]*parm[3]*(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 1 & d[[k-4]]$Phase == 2 & d[[k-4]]$trun == F)]^shape)) -
        sum(parm[k]*parm[3]*parm[4]*(d[[k-4]]$MBF[which(d[[k-4]]$Censor == 1 & d[[k-4]]$Phase == 3 & d[[k-4]]$trun == F)]^shape)) +
        sum(parm[k]*(d[[k-4]]$Lower[which(d[[k-4]]$Phase == 1 & d[[k-4]]$trun == T)]^shape) - 
              parm[k]*(d[[k-4]]$Upper[which(d[[k-4]]$Phase == 1 & d[[k-4]]$trun == T)]^shape)) +
        sum(parm[k]*parm[3]*(d[[k-4]]$Lower[which(d[[k-4]]$Phase == 2 & d[[k-4]]$trun == T)]^shape) - 
              parm[k]*parm[3]*(d[[k-4]]$Upper[which(d[[k-4]]$Phase == 2 & d[[k-4]]$trun == T)]^shape)) + 
        sum(parm[k]*parm[3]*parm[4]*(d[[k-4]]$Lower[which(d[[k-4]]$Phase == 3 & d[[k-4]]$trun == T)]^shape) - 
              parm[k]*parm[3]*parm[4]*(d[[k-4]]$Upper[which(d[[k-4]]$Phase == 3 & d[[k-4]]$trun == T)]^shape)) +
        dgamma(parm[k], parm[1], parm[2], log=T)
    }
    
    lp <- lp + 
      dgamma(parm[1], hyperA1, hyperA2, log=T) + 
      dgamma(parm[2],HyperG1, HyperG2, log=T) +
      dgamma(shape, shapePriorA, shapePriorB, log=T) +
      dgamma(parm[3],hyperT1A, hyperT1B, log=T) +
      dgamma(parm[4],hyperT2A, hyperT2B, log=T)
    
    return(lp)
  }
  
  logPostShape <- function(data, parms, shape, iter) {
    
    lpS <- 0
    
    for(w in 1:26){
      lpS <- lpS + logPostComp(data[[w]], parms[[w]][iter,], shape)
    }
    
    return(lpS)
  }
  
  # Save the starting data for imputation
  dataStart <- data
  
  draws <- list()
  
  for(i in 1:26){
    draws[[i]] <- matrix(0, nrow = samples, ncol = 12)
    
    draws[[i]][1,1] <- alphaStart
    draws[[i]][1,2] <- betaStart
    draws[[i]][1,3] <- theta1Start
    draws[[i]][1,4] <- theta2Start
  }

  draws[[27]] <- matrix(0, nrow = samples, ncol = 1)
  
  draws[[27]][1,1] <- shapeStart
  
  # phase obs counts 
  phase2Count <- rep(0, 26)
  phase3Count <- rep(0, 26)
  for(i in 1:26){
    for(j in 1:8){
      phase2Count[i] <- phase2Count[i] + length(which(data[[i]][[j]]$Phase == 2 & data[[i]][[j]]$Censor == 0))
      phase3Count[i] <- phase3Count[i] + length(which(data[[i]][[j]]$Phase == 3 & data[[i]][[j]]$Censor == 0))
    }
  }
  

  # counter for acceptance rate
  acca <- rep(0,26)
  accs <- 0

  for(i in 2:samples){

    for(j in 1:26){
      # update all 12 params
      
      for(k in 5:12){
        draws[[j]][i,k] <- rgamma(1,sum(data[[j]][[k-4]]$Censor == 0)+draws[[j]][i-1,1],
                             sum(data[[j]][[k-4]]$MBF[which(data[[j]][[k-4]]$Phase == 1)]^draws[[27]][i-1,1],
                                 (draws[[j]][i-1,3] * data[[j]][[k-4]]$MBF[which(data[[j]][[k-4]]$Phase == 2)]^draws[[27]][i-1,1]),
                                 (draws[[j]][i-1,3] * draws[[j]][i-1,4] * data[[j]][[k-4]]$MBF[which(data[[j]][[k-4]]$Phase == 3)]^draws[[27]][i-1,1]))
                             + draws[[j]][i-1,2])
      }
      
      draws[[j]][i,2] <- rgamma(1, 8*draws[[j]][i-1,1] + HyperG1,
                           sum(draws[[j]][i,5:12]) + HyperG2)
      
      # Phase sums
      phase2Sum <- 0
      phase3Sum3 <- 0
      phase3Sum4 <- 0
      for(k in 1:8){
        phase2Sum <- phase2Sum + draws[[j]][i,k+4] * sum(data[[j]][[k]]$MBF[which(data[[j]][[k]]$Phase == 2)]^draws[[27]][i-1,1])
        
        # for the first theta
        phase3Sum3 <- phase3Sum3 + draws[[j]][i,k+4] * draws[[j]][i-1,4] * sum(data[[j]][[k]]$MBF[which(data[[j]][[k]]$Phase == 3)]^draws[[27]][i-1,1])
        
        # for second theta
        phase3Sum4 <- phase3Sum4 + draws[[j]][i,k+4] * draws[[j]][i-1,3] * sum(data[[j]][[k]]$MBF[which(data[[j]][[k]]$Phase == 3)]^draws[[27]][i-1,1])
      }
      
      
      draws[[j]][i,3] <- rgamma(1, phase2Count[j] + phase3Count[j] + hyperT1A, phase2Sum + phase3Sum3 + hyperT1B)
      draws[[j]][i,4] <- rgamma(1, phase3Count[j] + hyperT2A, phase3Sum4 + hyperT2B)
      
      ################################
      # metropolis hastings
      ################################
      
      # Auto-adjusting Tuning Params
      if((i %% 500) == 0){
        if(acca[j]/i > .55){
          if(acca[j]/i > .7){
            tuningAV[j] <- tuningAV[j]*2
          } else{
            tuningAV[j] <- tuningAV[j]*1.5
          }
        } else if (acca[j]/i < .3) {
          if(acca[j]/i < .2){
            tuningAV[j] <- tuningAV[j]/2
          } else {
            tuningAV[j] <- tuningAV[j]/1.5
          }
        }
      }
      
      
      # Sample from alpha 
      draws[[j]][i,1] <- draws[[j]][i-1,1]
      astar <- rnorm(1, draws[[j]][i-1,1], tuningA)
      if (astar > 0) {
        lnew <- logPostComp(data[[j]], c(astar, draws[[j]][i,2:12]), draws[[27]][i-1,1])
        lold <- logPostComp(data[[j]], c(draws[[j]][i-1,1], draws[[j]][i,2:12]), draws[[27]][i-1,1])
        if(is.finite(lnew - lold)){
          if (lnew - lold > log(runif(1))) {
            draws[[j]][i,1] <- astar
            acca[j] <- acca[j] + 1
          }
        }
      }
    }
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      print(i)
      if(accs/i > .55){
        if(accs/i > .7){
          tuningS <- tuningS*2
        } else {
          tuningS <- tuningS*1.5
        }
      } else if (accs/i < .3) { 
        if (accs/i < .2){
          tuningS <- tuningS/2
        } else {
          tuningS <- tuningS/1.5
        }
      }
    }
    
    # try exponential
    if(EX == T){
      draws[[27]][i,1] <- 1
    } else {
      # update shape
      draws[[27]][i,1] <- draws[[27]][i-1,1]
      sstar <- rnorm(1, draws[[27]][i-1,1], tuningS)
      if (sstar > 0) {
        lnew <- logPostShape(data, draws, shape = sstar, iter = i)
        lold <- logPostShape(data, draws, shape = draws[[27]][i-1,1], iter = i)
        if(is.finite(lnew - lold)){
          if (lnew - lold > log(runif(1))) {
            draws[[27]][i,1] <- sstar
            accs <- accs + 1
          }
        }
      }
    }
  }
  
  devMat <- matrix(0, nrow = samples - 99, ncol = 26)
  
  # DIC
  if(samples > 200){
    d <- rep(0, samples-99)
    for(i in 100:samples){
      for(k in 1:26){
        for(j in 5:12){
          d[i-99] <- d[i-99] -
            2*(logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == F)], draws[[27]][i,1], draws[[k]][i,j]) +
                 logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == F)], draws[[27]][i,1], draws[[k]][i,j]*draws[[k]][i,3]) +
                 logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == F)], draws[[27]][i,1], draws[[k]][i,j]*draws[[k]][i,3]*draws[[k]][i,4]) -
                 sum(draws[[k]][i,j]*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == F)]^draws[[27]][i,1])) -
                 sum(draws[[k]][i,j]*draws[[k]][i,3]*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == F)]^draws[[27]][i,1])) -
                 sum(draws[[k]][i,j]*draws[[k]][i,3]*draws[[k]][i,4]*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == F)]^draws[[27]][i,1])) + 
                 sum(draws[[k]][i,j]*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1]) - 
                       draws[[k]][i,j]*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1])) +
                 sum(draws[[k]][i,j]*draws[[k]][i,3]*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1]) - 
                       draws[[k]][i,j]*draws[[k]][i,3]*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1])) + 
                 sum(draws[[k]][i,j]*draws[[k]][i,3]*draws[[k]][i,4]*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1]) - 
                       draws[[k]][i,j]*draws[[k]][i,3]*draws[[k]][i,4]*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == T)]^draws[[27]][i,1])))
        }
        if(k == 1){
          devMat[i-99,k] <- d[i-99]
        } else {
          devMat[i-99,k] <- d[i-99] - sum(devMat[i-99,1:k])
        }
      }
    }
    davg <- mean(d)

    dthetahat <- 0
    for(k in 1:26){
      for(j in 5:12){
        dthetahat <- dthetahat - 2*(logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == F)], mean(draws[[27]][100:samples,1]), mean(draws[[k]][100:samples,j])) +
                                      logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == F)], mean(draws[[27]][100:samples,1]), mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])) +
                                      logWeibullSumCPP(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 0 & data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == F)], mean(draws[[27]][100:samples,1]), mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*mean(draws[[k]][100:samples,4])) -
                                      sum(mean(draws[[k]][100:samples,j])*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == F)]^mean(draws[[27]][100:samples,1]))) -
                                      sum(mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == F)]^mean(draws[[27]][100:samples,1]))) -
                                      sum(mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*mean(draws[[k]][100:samples,4])*(data[[k]][[j-4]]$MBF[which(data[[k]][[j-4]]$Censor == 1 & data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == F)]^mean(draws[[27]][100:samples,1]))) +
                                      sum(mean(draws[[k]][100:samples,j])*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1])) - 
                                            mean(draws[[k]][100:samples,j])*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 1 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1]))) +
                                      sum(mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1])) - 
                                            mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 2 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1]))) + 
                                      sum(mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*mean(draws[[k]][100:samples,4])*(data[[k]][[j-4]]$Lower[which(data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1])) - 
                                            mean(draws[[k]][100:samples,j])*mean(draws[[k]][100:samples,3])*mean(draws[[k]][100:samples,4])*(data[[k]][[j-4]]$Upper[which(data[[k]][[j-4]]$Phase == 3 & data[[k]][[j-4]]$trun == T)]^mean(draws[[27]][100:samples,1]))))
      }
    }

    pd <- davg - dthetahat
    dic <- davg + pd
  } else {
    pd <- NULL
    dic <- NULL
    d <- NULL
  }

  return(list(draws = draws,
              acceptanceShape = accs/samples,
              acceptanceAlpha = acca/samples, 
              DIC = dic,
              PD = pd,
              Deviance = d,
              devMat = devMat))
}



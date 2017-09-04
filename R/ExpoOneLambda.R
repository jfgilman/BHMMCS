
library(coda)
library(LearnBayes)

# inverse CDF of weibull for generating data
eF.inv <- function(u,lam){
  return((-log(1-u)/lam))
}


ExpoOLMCMC <- function(data, samples, hyperB1 = .001, hyperB2 = .001,
                          hyperA1 = .001, hyperA2 = .001, alphaStart = 1,
                          hyperT1A = 3, hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                          betaStart = 1, theta1Start = 1, theta2Start = 1,
                          tuning = 1, trun = T){
  
  # Save the starting data for imputation
  dataStart <- data
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = 5)
  
  draws[1,1] <- alphaStart
  draws[1,2] <- betaStart
  draws[1,3] <- theta1Start
  draws[1,4] <- theta2Start
  
  # counter for acceptance rate
  acc <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    
    for(k in 1:length(d)){
      
      # 3rd line is subtracting sum of product of lamda and censored observations
      # this is the log likelihood of censored data
      lp <- lp + 
        sum(dexp(d[[k]][which(d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F),1], parm[5], log = T)) +
        sum(dexp(d[[k]][which(d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F),1], parm[5]*parm[3], log = T)) +
        sum(dexp(d[[k]][which(d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F),1], parm[5]*parm[3]*parm[4], log = T)) -
        sum(parm[5]*d[[k]][which(d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F),1]) -
        sum(parm[5]*parm[3]*d[[k]][which(d[[k]]$Censor == 1 & d[[k]]$Phase == 2 & d[[k]]$trun == F),1]) -
        sum(parm[5]*parm[3]*parm[4]*d[[k]][which(d[[k]]$Censor == 1 & d[[k]]$Phase == 3 & d[[k]]$trun == F),1]) + 
        sum(parm[k]*(d[[5]]$Lower[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)]) - 
              parm[k]*(d[[k]]$Upper[which(d[[k-5]]$Phase == 1 & d[[k]]$trun == T)])) +
        sum(parm[k]*parm[4]*(d[[k]]$Lower[which(d[[k-5]]$Phase == 2 & d[[k]]$trun == T)]) - 
              parm[k]*parm[4]*(d[[k]]$Upper[which(d[[k-5]]$Phase == 2 & d[[k]]$trun == T)])) + 
        sum(parm[k]*parm[4]*parm[5]*(d[[k]]$Lower[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)]) - 
              parm[k]*parm[4]*parm[5]*(d[[k]]$Upper[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)]))
    }
    
    lp <- lp + 
      dgamma(parm[1], hyperA1, hyperA2, log=T) + 
      dgamma(parm[2],hyperB1, hyperB2, log=T) + 
      dgamma(parm[3],hyperT1A, hyperT1B, log=T) + 
      dgamma(parm[4],hyperT2A, hyperT2B, log=T) +
      dgamma(parm[5], parm[1], parm[2], log=T)
    
    return(lp)
  }
  
  # phase obs counts 
  phase2Count <- 0
  phase3Count <- 0
  for(i in 1:length(data)){
    phase2Count <- phase2Count + length(which(data[[i]]$Phase == 2 & data[[i]]$Censor == 0))
    phase3Count <- phase3Count + length(which(data[[i]]$Phase == 3 & data[[i]]$Censor == 0))
  }
  
  totalObsNoCen <- sum(data[[1]]$Censor == 0,data[[2]]$Censor == 0,data[[3]]$Censor == 0,
                       data[[4]]$Censor == 0,data[[5]]$Censor == 0,data[[6]]$Censor == 0,
                       data[[7]]$Censor == 0,data[[8]]$Censor == 0)
  
  # MCMC draws
  for (i in 2:samples) {
    
    totalMBFP1 <- 0
    totalMBFP2 <- 0
    totalMBFP3 <- 0
    
    for(j in 1:length(data)){
      totalMBFP1 <- totalMBFP1 + sum(data[[j]]$MBF[which(data[[j]]$Phase == 1)])
      totalMBFP2 <- totalMBFP2 + sum(draws[i-1,3]*data[[j]]$MBF[which(data[[j]]$Phase == 2)])
      totalMBFP3 <- totalMBFP3 + sum(draws[i-1,3]*draws[i-1,4]*data[[j]]$MBF[which(data[[j]]$Phase == 3)])
    }
    
    draws[i,5] <- rgamma(1,totalObsNoCen+draws[i-1,1],
                           sum(totalMBFP1, totalMBFP2, totalMBFP3) + draws[i-1,2])

    draws[i,2] <- rgamma(1, draws[i-1,1] + hyperB1, draws[i,5] + hyperB2)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum3 <- 0
    phase3Sum4 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,5] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 2)])
      
      # for the first theta
      phase3Sum3 <- phase3Sum3 + draws[i,5] * draws[i-1,4] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
      
      # for second theta
      phase3Sum4 <- phase3Sum4 + draws[i,5] * draws[i-1,3] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
    }
    
    draws[i,3] <- rgamma(1, phase2Count + phase3Count + hyperT1A, phase2Sum + phase3Sum3 + hyperT1B)
    draws[i,4] <- rgamma(1, phase3Count + hyperT2A, phase3Sum4 + hyperT2B)
    
    ################################
    # metropolis hastings
    ################################
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      if(acc/i > .55){
        if(acc/i > .7){
          tuning <- tuning*2
        } else{
          tuning <- tuning*1.5
        }
      } else if (acc/i < .3) {
        if(acc/i < .2){
          tuning <- tuning/2
        } else {
          tuning <- tuning/1.5
        }
      }
    }
    
    #Sample from alpha 
    draws[i,1] <- draws[i-1,1]
    astar <- rnorm(1, draws[i-1,1], tuning)
    if (astar > 0) {
      lnew <- logPost(c(astar, draws[i,2:ncol(draws)]), data)
      lold <- logPost(c(draws[i-1,1], draws[i,2:ncol(draws)]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,1] <- astar
          acc <- acc + 1
        }
      }
    }
  }
  
  # need DIC
  
  # DIC
  if(samples > 200){
    d <- rep(0, samples-99)
    for(i in 100:samples){
      for(j  in 1:length(data)){
        d[i-99] <- d[i-99] - 
          2*(sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1)], draws[i,5], log = T)) +
               sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2)], draws[i,5]*draws[i,3], log = T)) +
               sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3)], draws[i,5]*draws[i,3]*draws[i,4], log = T)) -
               sum(draws[i,5]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1)])) -
               sum(draws[i,5]*draws[i,3]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2)])) -
               sum(draws[i,5]*draws[i,3]*draws[i,4]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3)])))
      }
    }
    
    davg <- mean(d)
    
    dthetahat <- 0
    for(j in 1:length(data)){
      dthetahat <- dthetahat - 2*(sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1)], mean(draws[100:samples,5]), log = T)) +
                                    sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2)], mean(draws[100:samples,5])*mean(draws[100:samples,3]), log = T)) +
                                    sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3)], mean(draws[100:samples,5])*mean(draws[100:samples,3])*mean(draws[100:samples,4]), log = T)) -
                                    sum(mean(draws[100:samples,5])*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1)])) -
                                    sum(mean(draws[100:samples,5])*mean(draws[100:samples,3])*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2)])) -
                                    sum(mean(draws[100:samples,5])*mean(draws[100:samples,3])*mean(draws[100:samples,4])*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3)])))
    }
    
    pd <- davg - dthetahat
    dic <- davg + pd
  } else {
    pd <- NULL
    dic <- NULL
  }
  
  
  return(list(draws = draws,
              acceptance = acc/samples,
              DIC = dic,
              PD = pd))
}


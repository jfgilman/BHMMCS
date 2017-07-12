library(Rcpp)

cppFunction('double logWeibullSumCPP(NumericVector x, double shape, double scale) {
            
            int nx = x.size();
            double d = 0.0;
            
            for(int i = 0; i < nx; i++) {
            d = d + log(shape) + log(scale) + (shape - 1)*log(x[i]) - scale*pow(x[i],shape);
            }
            
            return d;
            }')


# MCMC function for Weibull-Gamma Hierarchical Model
# Motivating example: JLTV a single component tested across multiple variants

WeibullPhaseMCMC <- function(data, samples = 5000, shapePriorA = .001,
                             shapePriorB = .001, HyperG1 = .001, HyperG2 = .001,
                             hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                             hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                             alphaStart = 1, betaStart = 1, theta1Start = 1,
                             theta2Start = 1, shapeStart = 1, tuningA = 1,
                             tuningS = 1){
  
  # Save the starting data for imputation
  dataStart <- data
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = length(data) + 5)
  
  draws[1,1] <- alphaStart
  draws[1,2] <- betaStart
  draws[1,3] <- shapeStart
  draws[1,4] <- theta1Start
  draws[1,5] <- theta2Start
  
  # counter for acceptance rate
  acca <- 0
  accs <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    
    for(k in 6:length(parm)){
      lp <- lp + 
        logWeibullSumCPP(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 1 & d[[k-5]]$trun == F)], parm[3], parm[k]) +
        logWeibullSumCPP(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 2 & d[[k-5]]$trun == F)], parm[3], parm[k]*parm[4]) +
        logWeibullSumCPP(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 3 & d[[k-5]]$trun == F)], parm[3], parm[k]*parm[4]*parm[5]) -
        sum(parm[k]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 1 & d[[k-5]]$trun == F)]^parm[3])) -
        sum(parm[k]*parm[4]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 2 & d[[k-5]]$trun == F)]^parm[3])) -
        sum(parm[k]*parm[4]*parm[5]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 3 & d[[k-5]]$trun == F)]^parm[3])) +
        sum(parm[k]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 1 & d[[k-5]]$trun == T)]^parm[3]) - 
              parm[k]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 1 & d[[k-5]]$trun == T)]^parm[3])) +
        sum(parm[k]*parm[4]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 2 & d[[k-5]]$trun == T)]^parm[3]) - 
              parm[k]*parm[4]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 2 & d[[k-5]]$trun == T)]^parm[3])) + 
        sum(parm[k]*parm[4]*parm[5]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 3 & d[[k-5]]$trun == T)]^parm[3]) - 
              parm[k]*parm[4]*parm[5]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 3 & d[[k-5]]$trun == T)]^parm[3])) +
        dgamma(parm[k], parm[1], parm[2], log=T)
    }
    
    lp <- lp + 
      dgamma(parm[1], hyperA1, hyperA2, log=T) + 
      dgamma(parm[2],HyperG1, HyperG2, log=T) +
      dgamma(parm[3], shapePriorA, shapePriorB, log=T) +
      dgamma(parm[4],hyperT1A, hyperT1B, log=T) +
      dgamma(parm[5],hyperT2A, hyperT2B, log=T)
    
    return(lp)
  }
  
  # phase obs counts 
  phase2Count <- 0
  phase3Count <- 0
  for(i in 1:length(data)){
    phase2Count <- phase2Count + length(which(data[[i]]$Phase == 2 & data[[i]]$Censor == 0))
    phase3Count <- phase3Count + length(which(data[[i]]$Phase == 3 & data[[i]]$Censor == 0))
  }
  
  # MCMC draws
  for (i in 2:samples) {
    
    for(j in 6:ncol(draws)){
      draws[i,j] <- rgamma(1,sum(data[[j-5]]$Censor == 0)+draws[i-1,1],
                           sum(data[[j-5]]$MBF[which(data[[j-5]]$Phase == 1)]^draws[i-1,3],
                               (draws[i-1,4] * data[[j-5]]$MBF[which(data[[j-5]]$Phase == 2)]^draws[i-1,3]),
                               (draws[i-1,4] * draws[i-1,5] * data[[j-5]]$MBF[which(data[[j-5]]$Phase == 3)]^draws[i-1,3]))
                           + draws[i-1,2])
    }
    
    draws[i,2] <- rgamma(1, length(data)*draws[i-1,1] + HyperG1,
                         sum(draws[i,6:ncol(draws)]) + HyperG2)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,k+5] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 2)]^draws[i-1,3])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + draws[i,k+5] * draws[i-1,5] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)]^draws[i-1,3])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + draws[i,k+5] * draws[i-1,4] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)]^draws[i-1,3])
    }
    
    draws[i,4] <- rgamma(1, phase2Count + phase3Count + hyperT1A, phase2Sum + phase3Sum4 + hyperT1B)
    draws[i,5] <- rgamma(1, phase3Count + hyperT2A, phase3Sum5 + hyperT2B)
    
    
    ################################
    # metropolis hastings
    ################################
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      if(acca/i > .55){
        if(acca/i > .7){
          tuningA <- tuningA*2
        } else{
          tuningA <- tuningA*1.5
        }
      } else if (acca/i < .3) {
        if(acca/i < .2){
          tuningA <- tuningA/2
        } else {
          tuningA <- tuningA/1.5
        }
      }
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
    
    # Sample from alpha 
    draws[i,1] <- draws[i-1,1]
    draws[i,3] <- draws[i-1,3]
    astar <- rnorm(1, draws[i-1,1], tuningA)
    if (astar > 0) {
      lnew <- logPost(c(astar, draws[i,2:ncol(draws)]), data)
      lold <- logPost(c(draws[i-1,1], draws[i,2:ncol(draws)]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,1] <- astar
          acca <- acca + 1
        }
      }
    }
    
    # Sample from shape
    sstar <- rnorm(1, draws[i-1,3], tuningS)
    if (sstar > 0) {
      lnew <- logPost(c( draws[i,1:2], sstar, draws[i,4:ncol(draws)]), data)
      lold <- logPost(c( draws[i,1:2], draws[i-1,3], draws[i,4:ncol(draws)]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,3] <- sstar
          accs <- accs + 1
        }
      }
    }
  }
  
  # DIC
  if(samples > 200){
    d <- rep(0, samples-99)
    for(i in 100:samples){
      for(j in 6:ncol(draws)){
        d[i-99] <- d[i-99] - 
          2*(logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)], draws[i,3], draws[i,j]) +
               logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)], draws[i,3], draws[i,j]*draws[i,4]) +
               logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)], draws[i,3], draws[i,j]*draws[i,4]*draws[i,5]) -
               sum(draws[i,j]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)]^draws[i,3])) -
               sum(draws[i,j]*draws[i,4]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)]^draws[i,3])) -
               sum(draws[i,j]*draws[i,4]*draws[i,5]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)]^draws[i,3])) +
               sum(draws[i,j]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]^draws[i,3]) - 
                     draws[i,j]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]^draws[i,3])) +
               sum(draws[i,j]*draws[i,4]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]^draws[i,3]) - 
                     draws[i,j]*draws[i,4]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]^draws[i,3])) + 
               sum(draws[i,j]*draws[i,4]*draws[i,5]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]^draws[i,3]) - 
                     draws[i,j]*draws[i,4]*draws[i,5]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]^draws[i,3])))
      }
    }
    
    davg <- mean(d)
    
    dthetahat <- 0
    for(j in 6:ncol(draws)){
      dthetahat <- dthetahat - 2*(logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)], mean(draws[100:samples,3]), mean(draws[100:samples,j])) +
                                    logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)], mean(draws[100:samples,3]), mean(draws[100:samples,j])*mean(draws[100:samples,4])) +
                                    logWeibullSumCPP(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)], mean(draws[100:samples,3]), mean(draws[100:samples,j])*mean(draws[100:samples,4])*mean(draws[100:samples,5])) -
                                    sum(mean(draws[100:samples,j])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)]^mean(draws[100:samples,3]))) -
                                    sum(mean(draws[100:samples,j])*mean(draws[100:samples,4])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)]^mean(draws[100:samples,3]))) -
                                    sum(mean(draws[100:samples,j])*mean(draws[100:samples,4])*mean(draws[100:samples,5])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)]^mean(draws[100:samples,3]))) +
                                    sum(mean(draws[100:samples,j])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3])) - 
                                          mean(draws[100:samples,j])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3]))) +
                                    sum(mean(draws[100:samples,j])*mean(draws[100:samples,4])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3])) - 
                                          mean(draws[100:samples,j])*mean(draws[100:samples,4])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3]))) + 
                                    sum(mean(draws[100:samples,j])*mean(draws[100:samples,4])*mean(draws[100:samples,5])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3])) - 
                                          mean(draws[100:samples,j])*mean(draws[100:samples,4])*mean(draws[100:samples,5])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]^mean(draws[100:samples,3]))))
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
              Deviance = d))
}
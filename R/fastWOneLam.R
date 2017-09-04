library(dplyr)
source("R/fastWLogLike.R")

# MCMC function for Weibull-Gamma Hierarchical Model
# Motivating example: JLTV a single component tested across multiple variants

fastWOneLam <- function(data, samples = 5000, shapePriorA = .001,
                      shapePriorB = .001, HyperG1 = .001, HyperG2 = .001,
                      hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                      hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                      alphaStart = 1, betaStart = 1, rho1Start = 1,
                      rho2Start = 1, shapeStart = 1, tuningA = 1,
                      tuningS = 1){
  
  # matrix for keeping MCMC draws for each parameter
  lam_draws <- rep(0, samples)
  lam_draws[1] <- 1
  
  alpha_draws <- rep(0, samples)
  alpha_draws[1] <- alphaStart
  
  beta_draws <- rep(0, samples)
  beta_draws[1] <- betaStart
  
  shape_draws <- rep(0, samples)
  shape_draws[1] <- shapeStart
  
  rho1_draws <- rep(0, samples)
  rho1_draws[1] <- rho1Start
  
  rho2_draws <- rep(0, samples)
  rho2_draws[1] <- rho2Start
  
  twoXlogLike <- rep(0, samples)
  
  # counter for acceptance rate
  acca <- 0
  accs <- 0
  
  # phase obs counts 
  phase2Count <- 0
  phase3Count <- 0
  for(i in 1:length(data)){
    phase2Count <- phase2Count + length(which(data[[i]]$Phase == 2 & data[[i]]$Censor == 0))
    phase3Count <- phase3Count + length(which(data[[i]]$Phase == 3 & data[[i]]$Censor == 0))
  }
  
  # MCMC draws
  for (i in 2:samples) {
    
    totalMBFP1 <- 0
    totalMBFP2 <- 0
    totalMBFP3 <- 0
    
    for(j in 1:length(data)){
      totalMBFP1 <- totalMBFP1 + sum(data[[j]]$MBF[which(data[[j]]$Phase == 1)]^shape_draws[i-1])
      totalMBFP2 <- totalMBFP2 + sum(rho1_draws[i-1]*(data[[j]]$MBF[which(data[[j]]$Phase == 2)]^shape_draws[i-1]))
      totalMBFP3 <- totalMBFP3 + sum(rho1_draws[i-1]*rho2_draws[i-1]*(data[[j]]$MBF[which(data[[j]]$Phase == 3)]^shape_draws[i-1]))
    }
    
    lam_draws[i] <- rgamma(1,totalObsNoCen + alpha_draws[i-1],
                           sum(totalMBFP1, totalMBFP2, totalMBFP3) + beta_draws[i-1])

    
    beta_draws[i] <- rgamma(1, length(data)*alpha_draws[i-1] + HyperG1,
                            sum(lam_draws[i,1:ncol(lam_draws)]) + HyperG2)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + lam_draws[i,k] * sum(data[[k]]$MBF[data[[k]]$Phase == 2]^shape_draws[i-1])
      
      # for the first rho
      phase3Sum4 <- phase3Sum4 + lam_draws[i,k] * rho2_draws[i-1] * sum(data[[k]]$MBF[data[[k]]$Phase == 3]^shape_draws[i-1])
      
      # for second rho
      phase3Sum5 <- phase3Sum5 + lam_draws[i,k] * rho1_draws[i-1] * sum(data[[k]]$MBF[data[[k]]$Phase == 3]^shape_draws[i-1])
    }
    
    rho1_draws[i] <- rgamma(1, phase2Count + phase3Count + hyperT1A, phase2Sum + phase3Sum4 + hyperT1B)
    rho2_draws[i] <- rgamma(1, phase3Count + hyperT2A, phase3Sum5 + hyperT2B)
    
    
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
    alpha_draws[i] <- alpha_draws[i-1]
    shape_draws[i] <- shape_draws[i-1]
    astar <- rnorm(1, alpha_draws[i-1], tuningA)
    if (astar > 0) {
      lnew <- logPost(d = data, lambdas = lam_draws[i], shape = shape_draws[i],
                      rho1 = rho1_draws[i], rho2 = rho2_draws[i],
                      alpha = astar, beta = beta_draws[i])
      lold <- logPost(d = data, lambdas = lam_draws[i], shape = shape_draws[i],
                      rho1 = rho1_draws[i], rho2 = rho2_draws[i],
                      alpha = alpha_draws[i-1], beta = beta_draws[i])
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          alpha_draws[i] <- astar
          acca <- acca + 1
        }
      }
    }
    
    # Sample from shape
    sstar <- rnorm(1, shape_draws[i-1], tuningS)
    if (sstar > 0) {
      lnew <- logPostOneLam(d = data, lambdas = lam_draws[i], shape = sstar,
                            rho1 = rho1_draws[i], rho2 = rho2_draws[i],
                            alpha = alpha_draws[i], beta = beta_draws[i])
      lold <- logPostOneLam(d = data, lambdas = lam_draws[i], shape = shape_draws[i-1],
                            rho1 = rho1_draws[i], rho2 = rho2_draws[i],
                            alpha = alpha_draws[i], beta = beta_draws[i])
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          shape_draws[i] <- sstar
          accs <- accs + 1
          
          twoXlogLike[i] <- (-2) * lnew
        } else {
          twoXlogLike[i] <- (-2) * lold
        }
      } else {
        twoXlogLike[i] <- (-2) * lold
      }
    }
  }
  
  # DIC
  if(samples > 200){
    
    davg <- mean(twoXlogLike[200:samples])
    
    dthetahat <- 0
    
    dthetahat <- logPost(d = data, lambdas = colMeans(lam_draws[200:samples]), 
                         shape = mean(shape_draws[200:samples]),
                         rho1 = mean(rho1_draws[200:samples]), 
                         rho2 = mean(rho2_draws[200:samples]),
                         alpha = mean(alpha_draws[200:samples]), 
                         beta = mean(beta_draws[200:samples]))
    
    pd <- davg + 2 * dthetahat
    dic <- davg + pd
  } else {
    pd <- NULL
    dic <- NULL
  }
  
  return(list(lam_draws = lam_draws,
              alpha_draws = alpha_draws,
              beta_draws = beta_draws,
              rho1_draws = rho1_draws,
              rho2_draws = rho2_draws,
              shape_draws = shape_draws,
              acceptanceShape = accs/samples,
              acceptanceAlpha = acca/samples, 
              DIC = dic,
              PD = pd))
  
}


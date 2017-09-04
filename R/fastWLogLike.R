
library(dplyr)
library(Rcpp)

cppFunction('double logWeibullSumCPP(NumericVector x, double shape, double scale) {

            int nx = x.size();
            double d = 0.0;
            
            for(int i = 0; i < nx; i++) {
            d = d + log(shape) + log(scale) + (shape - 1)*log(x[i]) - scale*pow(x[i],shape);
            }
            
            return d;
            }')

cppFunction('double logWeibullCenSumCPP(NumericVector x, double shape, double scale) {

            int nx = x.size();
            double d = 0.0;
            
            for(int i = 0; i < nx; i++) {
            d -= scale*pow(x[i],shape);
            }
            
            return d;
            }')

# log posterior function for MCMC draws
logPost <- function(d, lambdas, shape, rho1, rho2, alpha, beta, HyperG1 = .001,
                    HyperG2 = .001, hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                    hyperT1B = 3, hyperT2A = 3, hyperT2B = 3, shapePriorA = .001,
                    shapePriorB = .001) {
  
  lp <- 0
  
  for(k in 1:length(d)){
    
  lp <- lp + 
    logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambdas[k]) +
    logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F], shape, lambdas[k]*rho1) +
    logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F], shape, lambdas[k]*rho1*rho2) +
    logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambdas[k]) +
    logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambdas[k]*rho1) +
    logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambdas[k]*rho1*rho2) +
    sum(lambdas[k]*(d[[k]]$Lower[d[[k]]$Phase == 1 & d[[k]]$trun == T]^shape) - 
          lambdas[k]*(d[[k]]$Upper[d[[k]]$Phase == 1 & d[[k]]$trun == T]^shape)) +
    sum(lambdas[k]*rho1*(d[[k]]$Lower[d[[k]]$Phase == 2 & d[[k]]$trun == T]^shape) - 
          lambdas[k]*rho1*(d[[k]]$Upper[d[[k]]$Phase == 2 & d[[k]]$trun == T]^shape)) + 
    sum(lambdas[k]*rho1*rho2*(d[[k]]$Lower[d[[k]]$Phase == 3 & d[[k]]$trun == T]^shape) - 
          lambdas[k]*rho1*rho2*(d[[k]]$Upper[d[[k]]$Phase == 3 & d[[k]]$trun == T]^shape)) +
    dgamma(lambdas[k], alpha, beta, log=T)
  }
  
  lp <- lp + 
    dgamma(alpha, hyperA1, hyperA2, log=T) + 
    dgamma(beta,HyperG1, HyperG2, log=T) +
    dgamma(shape, shapePriorA, shapePriorB, log=T) +
    dgamma(rho1, hyperT1A, hyperT1B, log=T) +
    dgamma(rho2, hyperT2A, hyperT2B, log=T)
  
  return(lp)
}


# log posterior function for MCMC draws with one lambda
logPostOneLam <- function(d, lambda, shape, rho1, rho2, alpha, beta, HyperG1 = .001,
                    HyperG2 = .001, hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                    hyperT1B = 3, hyperT2A = 3, hyperT2B = 3, shapePriorA = .001,
                    shapePriorB = .001) {
  
  lp <- 0
  
  for(k in 1:length(d)){
    
    lp <- lp + 
      logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambda) +
      logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F], shape, lambda*rho1) +
      logWeibullSumCPP(d[[k]]$MBF[d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F], shape, lambda*rho1*rho2) +
      logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambda) +
      logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambda*rho1) +
      logWeibullCenSumCPP(d[[k]]$MBF[d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F], shape, lambda*rho1*rho2) +
      sum(lambda*(d[[k]]$Lower[d[[k]]$Phase == 1 & d[[k]]$trun == T]^shape) - 
            lambda*(d[[k]]$Upper[d[[k]]$Phase == 1 & d[[k]]$trun == T]^shape)) +
      sum(lambda*rho1*(d[[k]]$Lower[d[[k]]$Phase == 2 & d[[k]]$trun == T]^shape) - 
            lambda*rho1*(d[[k]]$Upper[d[[k]]$Phase == 2 & d[[k]]$trun == T]^shape)) + 
      sum(lambda*rho1*rho2*(d[[k]]$Lower[d[[k]]$Phase == 3 & d[[k]]$trun == T]^shape) - 
            lambda*rho1*rho2*(d[[k]]$Upper[d[[k]]$Phase == 3 & d[[k]]$trun == T]^shape)) +
      dgamma(lambda, alpha, beta, log=T)
  }
  
  lp <- lp + 
    dgamma(alpha, hyperA1, hyperA2, log=T) + 
    dgamma(beta,HyperG1, HyperG2, log=T) +
    dgamma(shape, shapePriorA, shapePriorB, log=T) +
    dgamma(rho1, hyperT1A, hyperT1B, log=T) +
    dgamma(rho2, hyperT2A, hyperT2B, log=T)
  
  return(lp)
}
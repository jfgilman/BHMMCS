
source("R/WeibullCDFfuns.R")

bayesChiSquare <- function(y, lambdas, betas, alpha = 0.95){
  
  n <- length(y)
  K <- n^0.4
  
  R_0 <- qchisq(alpha, K-1)
  
  R_B <- rep(0, length(lambdas))
  
  bins <- seq(from = 0, to = 1, by = 1/K)
  
  np <- n * (1/K)
  
  for(i in 1:length(lambdas)){

    F_lam_beta <- wF(y, lambdas[i], betas[i])
    for(j in 1:K){
      
      R_B[i] <- R_B[i] + ((sum(F_lam_beta >= bins[j] & F_lam_beta < bins[j + 1]) - np)^2) / np
    }
  }
  
  return( sum(R_B > R_0) / length(R_B))
  
}




test_y <- wF.inv(runif(100), 5, 3)
test_lam <- c(rep(0.2, 20), rep(0.04, 10), rep(0.005, 25))
test_beta <- c(rep(0.02, 20), rep(0.4, 10), rep(0.05, 25))

bayesChiSquare(test_y, test_lam, test_beta)



# CDF of weibull
wF <- function(x, scale, shape){
  return(1 - exp(-scale*(x)^shape))
}

# inverse CDF of weibull for generating data
wF.inv <- function(u,lam,beta){
  return((-log(1-u)/lam)^(1/beta))
}
source("WeibullCDFfuns.R")
source("WeibullPhaseMCMC.R")
source("fastWMcMC.R")


# simulation to test code
x <- vector("list", 8)

beta <- 0.7
lambda <- .2
theta1 <- .8
theta2 <- .9
for(i in 1:8){
  u <- runif(1000)
  x[[i]] <- data.frame(MBF = c(wF.inv(u[1:400],lambda,beta),
                               wF.inv(u[401:700],lambda*theta1,beta),
                               wF.inv(u[701:1000],lambda*theta1*theta2,beta)),
                       Censor = rep(0, 1000),
                       Phase = c(rep(1,400), rep(2,300), rep(3,300)),
                       trun = rep(F, 1000))
}


system.time({
  test1 <- WeibullPhaseMCMC(x, tuningA = 6, tuningS = .01,
                             samples = 5000, trun = F)
})

system.time({
  test2 <- fastWMcMC(x, tuningA = 6, tuningS = .01,
                             samples = 5000)
})

test1$acceptanceShape
test1$acceptanceAlpha

test2$acceptanceShape
test2$acceptanceAlpha

test1$DIC
test1$PD

test2$DIC
test2$PD

library(coda)

plot(as.mcmc(test1$draws[,1]))
plot(as.mcmc(test2$alpha_draws))

plot(as.mcmc(test1$draws[,2]))
plot(as.mcmc(test2$beta_draws))

plot(as.mcmc(test1$draws[,3]))
plot(as.mcmc(test2$shape_draws))

plot(as.mcmc(test1$draws[,4]))
plot(as.mcmc(test2$rho1_draws))

plot(as.mcmc(test1$draws[,5]))
plot(as.mcmc(test2$rho2_draws))

plot(as.mcmc(test1$draws[,6]))
plot(as.mcmc(test2$lam_draws[,1]))

plot(as.mcmc(test1$draws[,7]))
plot(as.mcmc(test2$lam_draws[,2]))











# CDF of weibull
wF <- function(x, scale, shape){
  return(1 - exp(-scale*(x)^shape))
}

# inverse CDF of weibull for generating data
wF.inv <- function(u,lam,beta){
  return((-log(1-u)/lam)^(1/beta))
}

# simulation test with no censoring or truncation 
shape <- 1
lambda <- rep(c(.2,.5,1,3), 8)
theta1 <- .8
theta2 <- .9

testData <- list()

for(i in 1:8){
  u <- runif(1000)
  testData[[i]] <- data.frame(MBF = c(wF.inv(u[1:400],lambda[i],shape),
                                           wF.inv(u[401:700],lambda[i]*theta1,shape),
                                           wF.inv(u[701:1000],lambda[i]*theta1*theta2,shape)),
                                   Censor = rep(0, 1000),
                                   Phase = c(rep(1,400), rep(2,300), rep(3,300)),
                                   trun = rep(F, 1000))
}

source("exponentialModel.R")

testResults <- ExpoPhaseMCMC(testData, samples = 5000)

testResults$acceptance

testResults$DIC
testResults$PD

library(coda)

plot(as.mcmc(testResults$draws[,1]))
plot(as.mcmc(testResults$draws[,2]))
plot(as.mcmc(testResults$draws[,3]))
plot(as.mcmc(testResults$draws[,4]))
plot(as.mcmc(testResults$draws[,5]))
plot(as.mcmc(testResults$draws[,6]))
plot(as.mcmc(testResults$draws[,7]))
plot(as.mcmc(testResults$draws[,8]))





# JLTV Data ---------------------------------------------------------------

JLTVData <- read.csv("data/JLTVData.csv")

JLTVData[, 1] <- as.numeric(JLTVData[, 1])
JLTVData$trun <- JLTVData$MBF == 0

JLTVData$Lower <- 0
JLTVData$Upper <- 0

subs <- list()

for(i in 1:26){
  subs[[i]] <- JLTVData[which(JLTVData[,9] == i),]
}

subList <- list()

for(i in 1:26){
  subList[[i]] <- list()
  for(j in 1:8){
    subList[[i]][[j]] <- subs[[i]][which(subs[[i]][,1]==j),c(6,8,3,10:12)]
  }
}


for(i in 1:length(subList)){
  for(j in 1:length(subList[[i]]))
    if(all(subList[[i]][[j]]$trun == F)){
    } else {
      for(k in 1:nrow(subList[[i]][[j]]))
        if(subList[[i]][[j]][k,4] == T){
          if(k == 1){
          } else {
            
            if(subList[[i]][[j]][k-1,2] == 1 & subList[[i]][[j]][k-1,3] != 3){
              subList[[i]][[j]][k,3] <- subList[[i]][[j]][k-1,3] + 1
            }
            
            if(subList[[i]][[j]][k+1,4] == F){
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+1,1]
            } else if(subList[[i]][[j]][k+2,4] == F){
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+2,1]
            } else if(subList[[i]][[j]][k+3,4] == F){
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+3,1]
            } else if(subList[[i]][[j]][k+4,4] == F){
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+4,1]
            } else if(subList[[i]][[j]][k+5,4] == F){
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+5,1]
            } else{
              subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+6,1]
            }
          }
        }
    }
}

testResults <- ExpoPhaseMCMC(subList[[2]], samples = 5000)

testResults$acceptance

testResults$DIC
testResults$PD

plot(as.mcmc(testResults$draws[,1]))
plot(as.mcmc(testResults$draws[,2]))
plot(as.mcmc(testResults$draws[,3]))
plot(as.mcmc(testResults$draws[,4]))
plot(as.mcmc(testResults$draws[,5]))
plot(as.mcmc(testResults$draws[,6]))
plot(as.mcmc(testResults$draws[,7]))
plot(as.mcmc(testResults$draws[,8]))

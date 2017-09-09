# Model Fitting
JLTVData <- read.csv("data/JLTVData.csv")

library(compiler)

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

iter <- 30000

# Exponential with some rate param from distribution and others common for each subsystem
source("R/ExpoOneLambda.R")
source("R/ExpoPhaseMCMC.R")
M2 <- cmpfun(ExpoOLMCMC)
M3 <- cmpfun(ExpoPhaseMCMC)

expoOutput <- list()

for(i in 1:26){
  if(i %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
    expoOutput[[i]] <- M2(subList[[i]], tuningA = 4, samples = iter)
  } else {
    expoOutput[[i]] <- M3(subList[[i]], tuningA = 4, samples = iter)
  }

  print(i)
}

save(expoOutput ,file="expoOutput.RData")

print("Expo Done")

# Weibull with some rate param from distribution and others common for each subsystem
source("R/FastWMcMC.R")
source("R/FastWOneLam.R")
M1 <- cmpfun(fastWMcMC)
M4 <- cmpfun(fastWOneLam)

weiOutput <- list()

for(i in 1:26){
  if(i %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
    weiOutput[[i]] <- M4(subList[[i]],
                       tuningA = 4,
                       tuningS = .02,
                       samples = iter)
  } else {
    weiOutput[[i]] <- M1(subList[[i]],
                       tuningA = 4,
                       tuningS = .02,
                       samples = iter)
  }
  print(i)
}

save(weiOutput ,file="weiOutput.RData")
print("Weibull Done")

# Weibull with a single shape parameter for all subsystems
# and all rate param from distribution

source("R/OneShapeMCMC.R")

M5 <- cmpfun(WeibullOneShapeMCMC)

output3 <- M5(subList, samples = iter, tuningS = .01, tuningA = 2)

save(output3 ,file="oneShapeOutput.RData")





# expoOutput[[1]]$acceptanceAlpha
# expoOutput[[2]]$acceptanceAlpha
# 
# expoOutput[[1]]$DIC
# expoOutput[[1]]$PD
# 
# expoOutput[[2]]$DIC
# expoOutput[[2]]$PD
# 
# library(coda)
# 
# plot(as.mcmc(expoOutput[[1]]$lam_draws))
# plot(as.mcmc(expoOutput[[2]]$lam_draws[,1]))
# 
# 
# 
# 
# weiOutput[[1]]$acceptanceShape
# weiOutput[[1]]$acceptanceAlpha
# 
# weiOutput[[2]]$acceptanceShape
# weiOutput[[2]]$acceptanceAlpha
# 
# weiOutput[[1]]$DIC
# weiOutput[[1]]$PD
# 
# weiOutput[[2]]$DIC
# weiOutput[[2]]$PD
# 
# library(coda)
# 
# plot(as.mcmc(weiOutput[[1]]$lam_draws))
# plot(as.mcmc(weiOutput[[2]]$lam_draws[,1]))




# Want to create test plan that controls for both consumer and producer risk
# Consumer risk: given test is passed, min prob expected failures is large
# P( E(fails in CMiles) > failLim | test passed) < alpha
#
# Producer risk: given test failed, min prob expected failures is small
# P( E(fails in PMiles) < failLim | test failed)
# 
# Cmiles < Pmiles
#
# vehNum: vehicle model number
# n: number of vehicles to test 
# allowedFailes: number of failures allowed to pass test


assuranceTest <- function(draws, vehNum = 1, n = 1, allowedFails = 0,
                          prodRisk = .05, conRisk = .1, failLim = 2,
                          Cmiles = 80, Pmiles = 140, startIndex = 1000,
                          mileStart = 500, increment = 25){
  
  miles <- mileStart
  
  running <- T
  while(running){
    
    if(miles %% 500 == 0){
      numeratorP <- 0
      denominatorP <- 0
      for(j in startIndex:nrow(draws[[1]])){
        rateParm <- 0
        for(i in 1:26){
          rateParm <- rateParm + 1/draws[[i]][j,3]*draws[[i]][j,4]*draws[[i]][j,vehNum + 4]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^draws[[27]][j,1])) * (rateParm*(miles^draws[[27]][j,1]))^k)/factorial(k)
        }
        if((Pmiles^draws[[27]][j,1])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      if(numeratorP/denominatorP > prodRisk){
        print("Producer risk criteria cannot be satisfied")
        print(numeratorP/denominatorP)
        return(NULL)
      }
    }
    
    print(miles)
    
    numeratorC <- 0
    denominatorC <- 0
    for(j in startIndex:nrow(draws[[1]])){
      rateParm <- 0
      for(i in 1:26){
        rateParm <- rateParm + 1/draws[[i]][j,3]*draws[[i]][j,4]*draws[[i]][j,vehNum + 4]
      }
      rateParm <- 1/rateParm
      amount <- 0
      for(k in 0:allowedFails){
        amount <- amount + (exp(-rateParm*(miles^draws[[27]][j,1])) * (rateParm*(miles^draws[[27]][j,1]))^k)/factorial(k)
      }
      if((Cmiles^draws[[27]][j,1])*rateParm > failLim){
        numeratorC <- numeratorC + amount
      }
      denominatorC <- denominatorC + amount
    }
    
    if(miles %% 500 == 0){
      print("Consumer risk alpha")
      print(numeratorC/denominatorC)
    }
    
    if(numeratorC/denominatorC < conRisk){
      numeratorP <- 0
      denominatorP <- 0
      for(j in startIndex:nrow(draws[[1]])){
        rateParm <- 0
        for(i in 1:26){
          rateParm <- rateParm + 1/draws[[i]][j,3]*draws[[i]][j,4]*draws[[i]][j,vehNum + 4]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^draws[[27]][j,1])) * (rateParm*(miles^draws[[27]][j,1]))^k)/factorial(k)
        }
        if((Pmiles^draws[[27]][j,1])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      
      if(numeratorP/denominatorP < prodRisk){
        running <- F
      } else {
        miles <- miles + increment
      }
    } else {
      miles <- miles + increment
    }
  }
  print("Final test miles")
  return(miles)
}


assuranceTest(output3$draws, vehNum = 1, failLim = 9, allowedFails = 1,
              startIndex = 1000, mileStart = 1000, increment = 100)



fails1 <- c()

for(j in 10000:nrow(output3$draws[[1]])){
  rateParm <- 0
  for(i in 1:26){
    rateParm <- rateParm + output3$draws[[i]][j,3]*output3$draws[[i]][j,4]*output3$draws[[i]][j,5]
  }
  fails1[j - 9999] <- (140^output3$draws[[27]][j,1])*rateParm 
}

hist(fails1)


fails2 <- c()

for(j in 10000:nrow(output3$draws[[1]])){
  rateParm <- 0
  for(i in 1:26){
    rateParm <- rateParm + output3$draws[[i]][j,3]*output3$draws[[i]][j,4]*output3$draws[[i]][j,5]
  }
  fails2[j - 9999] <- (80^output3$draws[[27]][j,1])*rateParm 
}

hist(fails2)

plot(density(fails2))
lines(density(fails1))



load("expoOutput.RData")

# i %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)

# phase 1 
phase1 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      gammab <- 1
       if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
         # gammab <- gamma((expoOutput[[j]]$shape_draws[i] + 1) / expoOutput[[j]]$shape_draws[i])
         holder[i-15000] <- expoOutput[[j]]$lam_draws[i]^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
         meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
       }else{
         holder[i-15000] <- expoOutput[[j]]$lam_draws[i,k]^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
         meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
       }

    }
    phase1[i-15000,k] <- 1/meanInvSum[i-15000]
  }
  
}

expo1 <- data.frame(rep("expo", 15000), rep("One", 15000), phase1[,1])

# for(i in 1:8){
#   print(mean(phase1[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.95)))
# }

# phase 2 
phase2 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      gammab <- 1
      if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
        holder[i-15000] <- (expoOutput[[j]]$lam_draws[i] * expoOutput[[j]]$rho1_draws[i])^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }else{
        holder[i-15000] <- (expoOutput[[j]]$lam_draws[i,k] * expoOutput[[j]]$rho1_draws[i])^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }
      
    }
    phase2[i-15000,k] <- 1/meanInvSum[i-15000]
  }
}

expo2 <- data.frame(rep("expo", 15000), rep("Two", 15000), phase2[,1])

# for(i in 1:8){
#   print(summary(phase1[,i]))
# }
# for(i in 1:8){
#   print(summary(phase2[,i]))
# }
# for(i in 1:8){
#   print(mean(phase2[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.95)))
# }

# phase 3 
phase3 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()

  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      gammab <- 1
      if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
        holder[i-15000] <- (expoOutput[[j]]$lam_draws[i] * expoOutput[[j]]$rho1_draws[i] * expoOutput[[j]]$rho2_draws[i])^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }else{
        holder[i-15000] <- (expoOutput[[j]]$lam_draws[i,k] * expoOutput[[j]]$rho1_draws[i] * expoOutput[[j]]$rho2_draws[i])^(-1/expoOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }
      
    }
    phase3[i-15000,k] <- 1/meanInvSum[i-15000]
  }
}

expo3 <- data.frame(rep("expo", 15000), rep("Three", 15000), phase3[,1])

# plot(density(phase3[,1]))
# for(i in 2:8){
#   lines(density(phase3[,i]))
# }
# 
# for(i in 1:8){
#   print(summary(phase3[,i]))
# }
# for(i in 1:8){
#   print(mean(phase3[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.95)))
# }

################################################################################

load("weiOutput.RData")

# i %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)

# phase 1 
phase1 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){

      if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- weiOutput[[j]]$lam_draws[i]^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }else{
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- weiOutput[[j]]$lam_draws[i,k]^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }
      
    }
    phase1[i-15000,k] <- 1/meanInvSum[i-15000]
  }
  
}

wei1<- data.frame(rep("Wei", 15000), rep("One", 15000), phase1[,1])

# for(i in 1:8){
#   print(mean(phase1[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.95)))
# }

# phase 2 
phase2 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      gammab <- 2
      if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- (weiOutput[[j]]$lam_draws[i] * weiOutput[[j]]$rho1_draws[i])^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }else{
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- (weiOutput[[j]]$lam_draws[i,k] * weiOutput[[j]]$rho1_draws[i])^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }
      
    }
    phase2[i-15000,k] <- 1/meanInvSum[i-15000]
  }
}

wei2<- data.frame(rep("Wei", 15000), rep("Two", 15000), phase2[,1])

# for(i in 1:8){
#   print(summary(phase2[,i]))
# }
# for(i in 1:8){
#   print(mean(phase2[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.95)))
# }

# phase 3 
phase3 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      gammab <- 2
      if(j %in% c(1, 4, 6, 7, 9, 12, 15, 16, 17, 18, 19, 21, 22, 24, 26)){
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- (weiOutput[[j]]$lam_draws[i] * weiOutput[[j]]$rho1_draws[i] * weiOutput[[j]]$rho2_draws[i])^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }else{
        gammab <- gamma((weiOutput[[j]]$shape_draws[i] + 1) / weiOutput[[j]]$shape_draws[i])
        holder[i-15000] <- (weiOutput[[j]]$lam_draws[i,k] * weiOutput[[j]]$rho1_draws[i] * weiOutput[[j]]$rho2_draws[i])^(-1/weiOutput[[j]]$shape_draws[i]) * gammab
        meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
      }
      
    }
    phase3[i-15000,k] <- 1/meanInvSum[i-15000]
  }
}

wei3<- data.frame(rep("Wei", 15000), rep("Three", 15000), phase3[,1])

# for(i in 1:8){
#   print(summary(phase3[,i]))
# }
# for(i in 1:8){
#   print(mean(phase3[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.95)))
# }

################################################################################

load("oneShapeOutput.RData")

# phase 1 
phase1 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    gammab <- gamma((output3$draws[[27]][i,1] + 1) / output3$draws[[27]][i,1])
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      holder[i-15000] <- output3$draws[[j]][i,k+4]^(-1/output3$draws[[27]][i,1]) * gammab
      meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
    }
    phase1[i-15000,k] <- 1/meanInvSum[i-15000]
  }
  
}

summary(output3$draws[[27]][,1])

wone1<- data.frame(rep("OneS", 15000), rep("One", 15000), phase1[,1])

for(i in 1:8){
  print(mean(phase1[,i]))
}
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase1[,i], c(.95)))
# }

# phase 2 
phase2 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    gammab <- gamma((output3$draws[[27]][i,1] + 1) / output3$draws[[27]][i,1])
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      holder[i-15000] <- (output3$draws[[j]][i,k+4]*output3$draws[[j]][i,3])^(-1/output3$draws[[27]][i,1]) * gammab
      meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
    }
    phase2[i-15000,k] <- 1/meanInvSum[i-15000]
  }
}

wone2<- data.frame(rep("OneS", 15000), rep("Two", 15000), phase2[,1])

# for(i in 1:8){
#   print(summary(phase2[,i]))
# }
# for(i in 1:8){
#   print(mean(phase2[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase2[,i], c(.95)))
# }

# phase 3 
phase3 <- matrix(0,nrow = 15000, ncol = 8)
for(k in 1:8){
  
  holder <- c()
  meanInvSum <- c()
  
  for(i in 15001:30000){
    gammab <- gamma((output3$draws[[27]][i,1] + 1) / output3$draws[[27]][i,1])
    meanInvSum[i-15000] <- 0
    for(j in 1:26){
      holder[i-15000] <- (output3$draws[[j]][i,k+4]*output3$draws[[j]][i,3]*output3$draws[[j]][i,4])^(-1/output3$draws[[27]][i,1]) * gammab
      meanInvSum[i-15000] <- meanInvSum[i-15000] + (1/holder[i-15000])
    }
    phase3[i-15000,k] <- 1/meanInvSum[i-15000]
  }
  
}

wone3 <- data.frame(rep("OneS", 15000), rep("Three", 15000), phase3[,1])

# for(i in 1:8){
#   print(summary(phase3[,i]))
# }
# for(i in 1:8){
#   print(mean(phase3[,i]))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.05)))
# }
# for(i in 1:8){
#   print(quantile(phase3[,i], c(.95)))
# }


colnames(expo1) <- c("Model", "Phase", "EMTF")
colnames(expo2) <- c("Model", "Phase", "EMTF")
colnames(expo3) <- c("Model", "Phase", "EMTF")
colnames(wei1) <- c("Model", "Phase", "EMTF")
colnames(wei2) <- c("Model", "Phase", "EMTF")
colnames(wei3) <- c("Model", "Phase", "EMTF")
colnames(wone1) <- c("Model", "Phase", "EMTF")
colnames(wone2) <- c("Model", "Phase", "EMTF")
colnames(wone3) <- c("Model", "Phase", "EMTF")
all.data <- rbind(expo1, expo2, expo3, wei1, wei2, wei3, wone1, wone2, wone3)

boxplots.triple = boxplot(all.data$EMTF ~ all.data$Phase + all.data$Model, at = c(1, 1.8, 2.6, 4, 4.8, 5.6, 7, 7.8, 8.6), xaxt='n', ylim = c(0, 800), col = c('white', 'white', 'gray'))
axis(side=1, at=c(1.8, 4.8, 7.8), labels=c('Exponential', 'Weibull', 'Weibull with One Shape'), line=0.5, lwd=0)
title('Comparing EMTF for Different Models')

rect(c(1.4, 4.4, 7.4), boxplots.triple$stats[2, c(2, 5, 8)], c(2.2, 5.2, 8.2), boxplots.triple$stats[4, c(2, 5, 8)], density=12, angle=45)
text(c(1, 1.8, 2.6, 4, 4.8, 5.6, 7, 7.8, 8.6), c(10, 10, 10, 10, 10, 10, 10, 10, 10), c('Phase 1', 'Phase 2', 'Phase 3', 'Phase 1', 'Phase 2', 'Phase 3', 'Phase 1', 'Phase 2', 'Phase 3'))





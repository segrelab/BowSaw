# July 26 
# Author: Demetrius DiMucci

# Sometimes there is a smaller network that gives you the same results
# This function takes the output of findNetwork() and evaluates all of the subnetworks.
# if it finds a subnetwork that gives the same confidence results it stops and returns that one to you

library(nsga2R)

curateNetwork <- function(N,target,Data,model,real){
  
  if(missing(real)){
    real = FALSE
  }
  if(real == TRUE){
    label <- model$y[target]
  } else {
    label <- model$predicted[target]
  }
  
  network <- N[[1]]
  
  matches <- seq(nrow(Data))
  for(j in 1:length(network)){
    matches <- intersect(matches,which(Data[,network[j]] == Data[target,network[j]]))
  }
  
  #print(length(matches))
  
  
  if(real == TRUE){
    initConf <- length(which(model$y[matches] == label))/length(matches)
  } else {
    initConf <- length(which(model$predicted[matches] == label))/length(matches)
  }
  #print(initConf)
  initMatches <- length(matches)
  
  # Knock out variables from the network - eliminate any variables that don't negatively affect
  # the results of the network. If all knocked out variables hurt network performance, stop and return
  # the remaining network.
  
  net = network
  
  prune = TRUE
  kill <- rep(0,length(net))
  
  while(prune == TRUE & length(net) > 1){
    
    
    
    kill <- rep(0,length(net))
    
    for(i in 1:length(net)){
      matches <- seq(nrow(Data))
      values <- Data[target,net]
      
      # Individually knockout variables and see how their absence affects network performance
      # If confidence is NOT negatively affected we know the variable is NON-Essential for making
      # the prediction and we designate it's index in the kill vector with a 1
      subnet <- net[-i]
      subVals <- values[-i]
      
      for(j in 1:length(subnet)){
        matches <- intersect(matches,which(Data[,subnet[j]] == subVals[j]))
      }
      
      
      if(real == FALSE){
        tempConf <- length(which(model$predicted[matches] == label))/length(matches)
        tempMatches <- length(matches)
      } else {
        tempConf <- length(which(model$y[matches] == label))/length(matches)
        tempMatches <- length(matches)
      }
      
      if(tempConf < initConf ){
        kill[i] <- 0
      } else {
        kill[i] <- 1
      }
    }
    if(sum(kill) > 0){
      net <- net[-which(kill == 1)]
    }
    
    if(sum(kill) < 1){
      prune = FALSE
    }
  }
  
  return(net)
}



# curateNetwork <- function(N,target,Data,RFclass,real){
#   model = RFclass
#   case = Data[target,]
#   network <- N[[1]]
#   
#   confidence <- 0
#   bestConf = max(N[[2]])
#   
#   count <- 0
#   while(confidence < bestConf){
#     count <- count + 1
#     tempNet <- combn(network,count)
#     
#     
#     tempConf <- vector()
#     tempMatches <- vector()
#     
#     for(j in 1:ncol(tempNet)){
#       matchSet <- seq(1:nrow(Data))
#       tempSet <- tempNet[,j]
#       for(k in 1:length(tempSet)){
#         matchSet <- intersect(matchSet,which(Data[,tempSet[k]] == Data[target,tempSet[k]]))
#       }
#       # If real == FALSE this will return a sub-network that predicts a greater fraction as the class
#       # If real == TRUE this will instead return a sub-network that has a greater fraction that is actually the predicted class
#       tempConf[j] <- length(which(model$predicted[matchSet] == model$predicted[target]))/length(matchSet)
#       if(real == TRUE){
#         tempConf[j] <- length(which(model$y[matchSet] == model$y[target]))/length(matchSet)
#       }
#       tempMatches[j] <- length(matchSet)
#     }
#     
#     confidence <- max(tempConf)
#   }
#   
#   effectiveNets <- which(tempConf >= bestConf)
#   results = list()
#   results[[1]]= t(tempNet[,effectiveNets])
#   results[[2]] = tempMatches[effectiveNets]
#   return(results)
#   
# }

#curateNetwork(N,target,Data,RFclass)

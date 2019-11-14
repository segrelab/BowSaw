# July 26 
# Author: Demetrius DiMucci

# This function takes the output of rankEdges() and finds an undirected small network that results in 100% confidence
# initRow indicates which edge of the ranked edge list to start from.

findSmallNetwork <- function(target,initRow,edges,Data,RFclass,real,categorical,minmatch){
  if(missing(minmatch)){
    minmatch = 10
  }
  if(missing(real)){
    real = FALSE
  }
  
  
  
  confidence = 0
  finalSet <- vector()
  keptEdges <- initRow + 1
  buddies <- vector()
  confCurve <- vector()
  model = RFclass
  DATA = Data
  binary <- length(unique(as.vector(DATA)))
  matchset <- seq(1:nrow(DATA))
  # Add edges to the network until either the forest is certain of its prediction or you run out of edges
  while(confidence < 1 & initRow < nrow(edges) ){
    initRow <- initRow + 1
    if(initRow > 1){
      featureSet <- edges[initRow,1:2]
    } else {
      featureSet <- edges[initRow,1:2]
    }
    
    if(categorical != TRUE){
      print('Data is not binary, recommended to use rankEdges_fuzzy()')
      if(initRow == 1){
        finalSet <- union(finalSet,featureSet)
        matchSet <- seq(1:nrow(DATA))
        for(k in 1:length(finalSet)){
          if(DATA[target,finalSet[k]] > medianDecision[finalSet[k]]){
            matchSet <- intersect(matchSet,which(DATA[,finalSet[k]] > medianDecision[finalSet[k]]))
          } else {
            matchSet <- intersect(matchSet,which(DATA[,finalSet[k]] <= medianDecision[finalSet[k]]))
          }
          
        }
        if(length(matchSet) >= minmatch){
          confidence <- length(which(model$predicted[matchSet] == model$predicted[target]))/length(matchSet)
          confCurve <- c(confCurve,confidence)
          buddies <- c(buddies,length(matchSet))
        }
        
      } else {
        matchSet <- seq(1:nrow(DATA))
        tempSet <- sort(union(finalSet,featureSet))
        for(k in 1:length(tempSet)){
          if(DATA[target,tempSet[k]] > medianDecision[tempSet[k]]){
            #print('more')
            matchSet <- intersect(matchSet,which(DATA[,tempSet[k]] > medianDecision[tempSet[k]]))
          } else {
            #print('less')
            matchSet <- intersect(matchSet,which(DATA[,tempSet[k]] <= medianDecision[tempSet[k]]))
          }
        }
        confidence2 <- length(which(model$predicted[matchSet] == model$predicted[target]))/length(matchSet)
        if(confidence2 > confidence & length(matchSet) >= minmatch){
          confidence = confidence2
          finalSet <- union(finalSet,featureSet)
          keptEdges <- c(keptEdges,initRow)
          buddies <- c(buddies,length(matchSet))
          confCurve <- c(confCurve,confidence)
        }
      }
    } else {
      if(initRow == 1){
        finalSet <- union(finalSet,featureSet)
        matchSet <- seq(1:nrow(DATA))
        for(k in 1:length(finalSet)){
          matchSet <- intersect(matchSet,which(DATA[,finalSet[k]] == DATA[target,finalSet[k]]))
        }
        confidence <- length(which(model$predicted[matchSet] == model$predicted[target]))/length(matchSet)
        if(real == TRUE){
          confidence <- length(which(model$y[matchSet] == model$y[target]))/length(matchSet)
        }
        buddies <- c(buddies,length(matchSet))
        confCurve <- c(confCurve,confidence)
      } else {
        matchSet <- seq(1:nrow(DATA))
        tempSet <- union(finalSet,featureSet)
        for(k in 1:length(tempSet)){
          matchSet <- intersect(matchSet,which(DATA[,tempSet[k]] == DATA[target,tempSet[k]]))
        }
        confidence2 <- length(which(model$predicted[matchSet] == model$predicted[target]))/length(matchSet)
        if(real == TRUE){
          confidence2 <- length(which(model$y[matchSet] == model$y[target]))/length(matchSet)
        }
        
        
        if(confidence2 > confidence & length(matchSet) >= minmatch){
          confidence = confidence2
          finalSet <- union(finalSet,featureSet)
          keptEdges <- c(keptEdges,initRow)
          buddies <- c(buddies,length(matchSet))
          confCurve <- c(confCurve,confidence)
        }
        
      }
    }
  }
  results <- list()
  results[[1]] = finalSet
  results[[2]] = confCurve
  results[[3]] = keptEdges
  results[[4]] = buddies
  return(results)
}

#N = findSmallNetwork(target,6,newEDGES,X,RFclass)



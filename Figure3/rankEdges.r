# July 26 
# Author: Demetrius DiMucci

# This function takes the output of countNetwork() and generates a ranked list of edges
# Default is the edges are sorted by their frequency. Additional metrics are provided for potential optimization later

rankEdges <- function(target,Combinations,DATA,RFclass,threshold,categorical){
  model = RFclass
  
  if(missing(categorical)){
    categorical = FALSE
  }
  
  #binary <- length(unique(as.vector(DATA)))
  if(missing(threshold)){
    threshold = 3
  }
  
  if(!(categorical)){
    print('The features are not categorical, will not be able to calculate confidences. Use functions X then Y')
  } else {
    count <- 0
    # To save computational time we consider only those edges that appeared more than the threshold setting
    check <- which(Combinations > threshold, arr.ind = TRUE)
    EDGES <- matrix(0,ncol=9,nrow=nrow(check))
    colnames(EDGES) <- c('Proximate','Prior','Pval','Matches','Confidence','Association','Frequency','val1','val2')
    for(k in 1:nrow(check)){
      count = k
      EDGES[count,1] = check[k,1]
      EDGES[count,2] = check[k,2]
      
      i = check[k,1]
      j = check[k,2]
      
      # Calculate the probability of seeing a given j -> i edge more often than we did
      DrawnWhiteBalls <- Combinations[i,j]
      totalBalls <- sum(Combinations)
      whiteBalls <- sum(Combinations[i,])
      blackBalls <- totalBalls - whiteBalls
      drawnBalls <- sum(Combinations[,j])
      EDGES[count,3] = phyper(DrawnWhiteBalls-1,whiteBalls,blackBalls,drawnBalls,lower.tail=FALSE)
      
      # Identify how many other samples have the same values as the target sample at this edge
      # This currently includes the target sample
      m <- which(DATA[,i] == DATA[target,i] & DATA[,j] == DATA[target,j])
      b <- which(DATA[,i] == DATA[target,i] & DATA[,j] == DATA[target,j] & model$y == 1)
      EDGES[count,4] <- length(m)
      
      # If we were to split the samples based only on this edge 
      # how confident is the forest that the sample is the final predicted class?
      EDGES[count,5] <- length(which(model$predicted[m] == model$predicted[target]))/length(m)
      # How strongly is the edge associated with the observed label?
      EDGES[count,6] <- length(which(model$y[m] == model$y[target]))/length(m)
      EDGES[count,7] <- DrawnWhiteBalls # How many times this particular edge occurred
      EDGES[count,8] <-  DATA[target,i]
      EDGES[count,9] <-  DATA[target,j]
    }
    
    newEDGES <- EDGES[rev(order(EDGES[,7])),]
    return(newEDGES)
  }
}

#newEDGES <- rankEdges(S,Data,RFclass)

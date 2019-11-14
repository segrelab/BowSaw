prioritize <- function(cand1,cand2,x,mod,lab){
  labs = mod$y
  cand3  = cand1
  bestOrd <- rep(0,nrow(cand3))
  bestOrd[1] = 1
  accounted <- bestOrd*0
  
  
  RuleDist <- matrix(0,ncol=3,nrow=nrow(cand3))
  rule1 = which(cand3[1,]>0)
  vals1 = cand3[1,rule1]
  
  rule1Matches <- seq(nrow(x))
  for(i in 1:length(rule1)){
    rule1Matches <- intersect(rule1Matches,which(x[,rule1[i]] == vals1[i]))
  }
  RuleDist[1,3] = length(rule1Matches)
  allMatches <- rule1Matches
  
  added <- vector()
  added[1] <- length(rule1Matches)
  accounted[1] = length(rule1Matches)/length(which(labs == lab))
  for(k in 2:nrow(cand3)){
    if(sum(added) < length(which(labs == lab))){
      
      newAdditions <- vector()
      
      for(i in 1:nrow(cand3)){
        
        candRule <- which(cand3[i,]>0)
        candVals <- cand3[i,candRule]
        
        candMatches <- seq(nrow(x))
        for(j in 1:length(candRule)){
          candMatches <- intersect(candMatches,which(x[,candRule[j]] == candVals[j]))
        }
        
        newAdditions[i] <- (length(candMatches) - length(intersect(candMatches,allMatches)))
        
      }
      bestOrd[k] <- which.max(newAdditions)
      added[k] <- newAdditions[bestOrd[k]]
      candRule <- which(cand3[bestOrd[k],]>0)
      candVals <- cand3[bestOrd[k],candRule]
      
      candMatches <- seq(nrow(x))
      for(j in 1:length(candRule)){
        candMatches <- intersect(candMatches,which(x[,candRule[j]] == candVals[j]))
      }
      allMatches <- union(allMatches,candMatches)
      accounted[k] <- length(allMatches)/length(which(labs == lab))
    }
  }
  
  
  alls = seq(nrow(cand1))
  z <- max(which(bestOrd > 0))
  zz <- alls[bestOrd[1:z]]
  newAlls <- c(alls[zz],alls[-zz])
  
  newCand <- cand1[newAlls,]
  newCand2 <- cand2[newAlls,]
  
  results <- list()
  results[[1]] = newCand
  results[[2]] = newCand2
  results[[3]] = bestOrd
  results[[4]] = accounted
  results[[5]] = added
  return(results)
}

######################################################33
unchanged <- function(cand1,cand2,x,mod,lab){
  labs = mod$y
  cand3  = cand1
  bestOrd <- seq(nrow(cand3))
  accounted <- bestOrd*0
  
  
  RuleDist <- matrix(0,ncol=3,nrow=nrow(cand3))
  rule1 = which(cand3[1,]>0)
  vals1 = cand3[1,rule1]
  
  rule1Matches <- seq(nrow(x))
  for(i in 1:length(rule1)){
    rule1Matches <- intersect(rule1Matches,which(x[,rule1[i]] == vals1[i]))
  }
  RuleDist[1,3] = length(rule1Matches)
  allMatches <- rule1Matches
  
  added <- vector()
  added[1] <- length(rule1Matches)
  accounted[1] = length(rule1Matches)/length(which(labs == lab))
  for(k in 2:nrow(cand3)){
    candRule <- which(cand3[k,]>0)
    candVals <- cand3[k,candRule]
    
    candMatches <- seq(nrow(x))
    for(j in 1:length(candRule)){
      candMatches <- intersect(candMatches,which(x[,candRule[j]] == candVals[j]))
    }
    added[k] <- (length(candMatches) - length(intersect(candMatches,allMatches)))
    allMatches <- union(allMatches,candMatches)
  }
  
  
  
  results <- list()
  results[[1]] = cand1
  results[[2]] = cand2
  results[[3]] = bestOrd
  results[[4]] = cumsum(added/length(which(labs==lab)))
  results[[5]] = added
  return(results)
}
###########################################################################################

# This function is the same as prioirtize except that it allows for imperfect matches (association < 1)
prioritize_imperfect <- function(cand1,cand2,x,mod,lab){
  labs = mod$y
  total <- length(which(labs == lab))
  cand3  = cand1
  bestOrd = 1
  accounted <- bestOrd*0
  totAdded <- vector()
  labeled <- vector()
  
  
  RuleDist <- matrix(0,ncol=3,nrow=nrow(cand3))
  rule1 = which(cand3[1,]>0)
  vals1 = cand3[1,rule1]
  
  rule1Matches <- seq(nrow(x))
  for(i in 1:length(rule1)){
    rule1Matches <- intersect(rule1Matches,which(x[,rule1[i]] == vals1[i]))
  }
  RuleDist[1,3] = length(rule1Matches)
  allMatches <- rule1Matches
  realMatches <- allMatches[which(labs[allMatches]== lab)]
  
  added <- vector()
  added[1] <- length(which(labs[realMatches]==lab))
  accounted[1] = length(realMatches)/length(which(labs == lab))
  totAdded[1] <- length(realMatches)
  labeled[1] <- added[1]
  #unlabeled[1] <- 0
  ##############################
  cycle = 1
  while(sum(added) < total & cycle < 200){
    cycle <- cycle + 1
    newAdds <- vector()
    for(k in 1:nrow(cand3)){
      candRule <- which(cand3[k,]>0)
      candVals <- cand3[k,candRule]
      
      candMatch <- seq(nrow(x))
      for(j in 1:length(candRule)){
        candMatch <- intersect(candMatch,which(x[,candRule[j]] == candVals[j]))
      }
      
      goodCand <- candMatch[which(labs[candMatch] == lab)]
      tempMatch <- union(goodCand,realMatches)
      newAdds[k] <- length(tempMatch) - length(realMatches)
    }
    best <- which.max(newAdds)
    bestOrd <- c(bestOrd,best)
    candRule <- which(cand3[best,]>0)
    candVals <- cand3[best,candRule]
    
    candMatch <- seq(nrow(x))
    for(j in 1:length(candRule)){
      candMatch <- intersect(candMatch,which(x[,candRule[j]] == candVals[j]))
    }
    goodCand <- candMatch[which(labs[candMatch] == lab)]
    badCand <- candMatch[which(labs[candMatch] != lab)]
    
    tempMatch <- union(goodCand,realMatches)
    alltemp <- union(candMatch,allMatches)
    added <- c(added,length(tempMatch) - length(realMatches))
    totAdded <- c(totAdded,length(alltemp) - length(allMatches))
    labeled <- c(labeled,length(goodCand))
    #unlabeled <- c(labeled,length(badCand))
    
    allMatches <- union(allMatches,candMatch)
    realMatches <- union(realMatches,goodCand)
  }
  
  
  #####################
  newCand <- cand1[bestOrd,]
  newCand2 <- cand2[bestOrd,]
  
  results <- list()
  results[[1]] = newCand
  results[[2]] = newCand2
  results[[3]] = bestOrd
  results[[4]] = totAdded
  results[[5]] = added
  results[[6]] = labeled
  #results[[7]] = unlabeled
  return(results)
}

# This function is the same as prioirtize except that it allows for imperfect matches (association < 1)
prioritize_imperfect <- function(cand1,cand2,x,mod,lab){
  labs = mod$y
  total <- length(which(labs == lab))
  cand3  = cand1
  bestOrd = 1
  accounted <- bestOrd*0
  totAdded <- vector()
  labeled <- vector()
  
  
  RuleDist <- matrix(0,ncol=3,nrow=nrow(cand3))
  rule1 = which(cand3[1,]>0)
  vals1 = cand3[1,rule1]
  
  rule1Matches <- seq(nrow(x))
  for(i in 1:length(rule1)){
    rule1Matches <- intersect(rule1Matches,which(x[,rule1[i]] == vals1[i]))
  }
  RuleDist[1,3] = length(rule1Matches)
  allMatches <- rule1Matches
  realMatches <- allMatches[which(labs[allMatches]== lab)]
  
  added <- vector()
  added[1] <- length(which(labs[realMatches]==lab))
  accounted[1] = length(realMatches)/length(which(labs == lab))
  totAdded[1] <- length(realMatches)
  labeled[1] <- added[1]
  #unlabeled[1] <- 0
  ##############################
  cycle = 1
  while(sum(added) < total & cycle < 2000){
    cycle <- cycle + 1
    newAdds <- vector()
    for(k in 1:nrow(cand3)){
      candRule <- which(cand3[k,]>0)
      candVals <- cand3[k,candRule]
      
      candMatch <- seq(nrow(x))
      for(j in 1:length(candRule)){
        candMatch <- intersect(candMatch,which(x[,candRule[j]] == candVals[j]))
      }
      
      goodCand <- candMatch[which(labs[candMatch] == lab)]
      tempMatch <- union(goodCand,realMatches)
      newAdds[k] <- length(tempMatch) - length(realMatches)
    }
    best <- which.max(newAdds)
    bestOrd <- c(bestOrd,best)
    candRule <- which(cand3[best,]>0)
    candVals <- cand3[best,candRule]
    
    candMatch <- seq(nrow(x))
    for(j in 1:length(candRule)){
      candMatch <- intersect(candMatch,which(x[,candRule[j]] == candVals[j]))
    }
    goodCand <- candMatch[which(labs[candMatch] == lab)]
    badCand <- candMatch[which(labs[candMatch] != lab)]
    
    tempMatch <- union(goodCand,realMatches)
    alltemp <- union(candMatch,allMatches)
    added <- c(added,length(tempMatch) - length(realMatches))
    totAdded <- c(totAdded,length(alltemp) - length(allMatches))
    labeled <- c(labeled,length(goodCand))
    #unlabeled <- c(labeled,length(badCand))
    
    allMatches <- union(allMatches,candMatch)
    realMatches <- union(realMatches,goodCand)
  }
  
  
  #####################
  newCand <- cand1[bestOrd,]
  newCand2 <- cand2[bestOrd,]
  
  results <- list()
  results[[1]] = newCand
  results[[2]] = newCand2
  results[[3]] = bestOrd
  results[[4]] = totAdded
  results[[5]] = added
  results[[6]] = labeled
  #results[[7]] = unlabeled
  return(results)
}




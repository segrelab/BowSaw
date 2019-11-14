# 20 November 18
# Author: Demetrius DiMucci

# This function will:
# Take the output of BowSaw() and rank each of the rules by their appropriateness as a true rule.

# Helper function
# Starting at the rule with the most hits we will progress down the list of candidate rules
# and evaluate if the rule is a subset of each rule - if it is then we will ignore the larger rule
# as it will have fewer matches. If not we will not mark it. We will do this recursively until
# no more rules can be marked ignore. These will be the candidate rules.
isSubset <- function(rule1,rule2,rul1,rul2){
  
  
  if(length(rule1) > 0 & length(rule2) > 0){
    subset = length(which(rule1 %in% rule2))/length(rule1)
    
    if(subset == 1){
      check <- which(rule2 %in% rule1)
      dif = 0
      for(i in 1:length(check)){
        if(rul2[check[i]] != rul1[i]){
          dif = dif + 1
        }
      }
      if(dif != 0){
        subset = 0
      }
    }
  } else {
    subset = 0
  }
  
  return(subset)
}

#####
# evaluateBowSaw will return:
# [[1]] - the matrix of Rules encountered by BowSaw 
# [[2]] - how many samples contain that rule and if the rule contains a complete rule discovered by bowsaw
# The list is returned in descending order of number of matches - naive assumption is that more matches is a better rule
evaluateBowSaw <- function(RULES,x,mod){
  labs = mod$y
  # Record how long each rule is - we start from the smallest rules and move up if necessary
  # Record how many samples match to each rule - when everything has full penetrance this is the only metric
  # we have to evaluate goodness
  rights <- vector()
  matches <- vector()
  for(i in 1:nrow(RULES)){
    rule <- which(RULES[i,] != 0)
    vals <- RULES[i,rule]
    
    if(length(rule) > 1){
      
      netMatches <- seq(nrow(x))
      for(j in 1:length(rule)){
        netMatches <- intersect(netMatches,which(x[,rule[j]] == vals[j]))
      }
      matches[i] <- length(netMatches)
      rights[i] <- length(which(labs[netMatches] == 1))
    } else {
      netMatches <- which(x[,rule] ==vals)
      matches[i] <- length(netMatches)
      rights[i] <- length(which(labs[netMatches] == 1))
    }
  }
  
  ignore <- rep(0,nrow(RULES))
  candidates <- cbind(matches,rights,ignore)
  
  ord <- rev(order(matches))
  RULES2 <- RULES[ord,]
  candidates <- candidates[ord,]
  lengths <- vector()
  
  for(i in 1:nrow(candidates)){
    rule <- which(RULES2[i,] != 0)
    lengths[i] <- length(rule)
    rul1 <- RULES2[i,rule]
    if(candidates[i,2] == 0){
      for(j in 1:nrow(candidates)){
        if(i != j){
          rule2 <- which(RULES2[j,] != 0)
          rul2 <- RULES2[j,rule2]
          
          subset <- isSubset(rule,rule2,rul1,rul2)
          #print(subset)

          if(subset == 1){
            candidates[j,2] = 1
          }
        }
      }
    }
  }
  candidates <- cbind(candidates,lengths)
  colnames(candidates) <- c('Matches','Labeled 1','Contains Subset','Rule Length')
  results <- list()
  results[[1]] <- RULES2
  results[[2]] <- candidates
  return(results)
}



####
# This function will evaluate each candidate rule and ask what happens if we prune the rule
evaluateCandidates <- function(candidates,RULES2,mod){
  labs <- mod$y
  targets <- which(candidates[,2] == 0)
  pruned <- list()
  change <- vector()
  ruleLen <- vector()
  accs <- vector()
  newAccs <- vector()
  newMatch <- vector()
  for(i in 1:length(targets)){
    rule = which(RULES2[targets[i],] != 0)
    vals <- RULES2[targets[i],rule]
    if(length(rule) > 1){
      netMatches <- seq(nrow(x))
      for(j in 1:length(rule)){
        netMatches <- intersect(netMatches,which(x[,rule[j]] == vals[j]))
      }
      match = netMatches
    } else {
      match <- which(x[,rule] == vals)
    }
    
    initMatch <- length(match)
    initAcc <- length(which(labs[match] == 1))/length(match)
    accs[i] = initAcc
    
    tempRules <- list()
    tempAcc <- vector()
    tempMatch <- vector()
    for(j in 1:length(rule)){
      if(length(rule) != 2){
        vals <- (RULES2[i,rule[-j]])
        tempRule <- rule[-j]
        netMatches <- seq(nrow(x))
        for(k in 1:length(tempRule)){
          netMatches <- intersect(netMatches,which(x[,tempRule[k]] == vals[k]))
        }
        match = netMatches
        tempAcc[j] <- length(which(labs[match] == 1))/length(match)
        tempMatch[j] <- length(match)
        tempRules[[j]] <- rule[-j]
      }
    }
    
    if(length(which(is.na(tempAcc))) > 0){
      tempAcc[which(is.na(tempAcc))] = 0
    }
    
    better <- which(tempAcc >= initAcc)
    
    if(length(better) > 0){
      change[i] = 'TRUE'
      best <- which.max(tempAcc)
      newAccs[i] <- tempAcc[best]
      newMatch[i] <- tempMatch[best]
      pruned[[i]] <- rule[-best]
      ruleLen[i] <- length(pruned[[i]])
    } else {
      change[i] = FALSE
      newAccs[i] <- initAcc
      newMatch[i] <- 0
      pruned[[i]] <- rule
      ruleLen[i] <- length(pruned[[i]])
    }
  }
  results = list()
  results[[1]] = pruned
  results[[2]] = cbind(targets,ruleLen,candidates[targets,],accs,newMatch,newAccs,change)
  results[[3]] = RULES2
  results[[4]] = RULES2[targets,]

  return(results)
}


#### fraction represented
# This function will go down the list of candidate rules and keep a running tally of how many
# samples with the observed label are accounted for. If all of the samples are accounted for by all the candidate
# rules this will suggest that pruning is unnecessary.

labelAccounting <- function(cand1,x,labs){
  results <- list()
  
  keptRules <- rep(0,nrow(cand1))
  allSamps <- vector()
  sampsFound <- rep(0,nrow(cand1))
  labsFound <- rep(0,nrow(cand1))
  targetMat <- matrix(0,ncol=nrow(cand1),nrow=nrow(x))
  addedValue <- rep(0,nrow(cand1))
  
  for(i in 1:nrow(cand1)){
    rule <- which(cand1[i,]!=0)
    vals <- cand1[i,rule]
    
    matches <- seq(nrow(x))
    for(j in 1:length(rule)){
      matches <- intersect(matches,which(x[,rule[j]] == vals[j]))
    }
    
    if(length(matches) > 0){
      good <- which(labs[matches] == 1 )
      targetMat[matches[good],i] = 1
    }
    
    alreadyFound <- which(matches %in% allSamps)
    if(length(alreadyFound) != length(matches)){
      keptRules[i] = 1
      addedValue[i] <- length(matches) - length(alreadyFound)
    }
    
    allSamps <- union(allSamps,matches)
    sampsFound[i] <- length(allSamps)
    labsFound[i] <- length(which(labs[allSamps] == 1))
    
    
    # if(labsFound[i]/length(which(labs ==1)) == 1){
    #   break()
    # }
  }
  results[[1]] <- allSamps
  results[[2]] <- sampsFound
  results[[3]] <- labsFound
  results[[4]] <- labsFound/length(which(labs ==1))
  results[[5]] <- cand1[which(keptRules == 1),]
  results[[6]] <- keptRules
  results[[7]] <- targetMat
  results[[8]] <- addedValue
  return(results)
}


##### makePlot() is a funciton that will take all candidate rules that combined account for more than a 
# user defined threshold of the labeled data and draw a heatmap of the jaccard similarities between them
# The intention is to use this visualization as a means to estimate how many underlying rules exist in the
# data
makePlot <- function(cand1,cand2,targs,threshold,counts){
  cand3 <- cand1[which(counts <= threshold & counts > 0),]
  JD <- matrix(0,ncol=nrow(cand3),nrow=nrow(cand3))
  
  nameList <- vector()
  for(i in 1:nrow(cand3)){
    realName <- vector()
    temp <- paste0(as.character(which(cand3[i,]>0)))
    for(k in 1:length(temp)){
      realName <- paste0(realName,temp[k],',')
    }
    nameList[i] <- realName
  }
  colnames(JD) = nameList
  rownames(JD) = nameList
  for(i in 1:nrow(JD)){
    for(j in 1:nrow(JD)){
      rule1 <- which(targs[,i]>0)
      rule2 <- which(targs[,j] > 0)
      
      dist <- length(intersect(rule1,rule2))/length(union(rule1,rule2))
      JD[i,j] = dist
      JD[j,i] = dist
    }
  }
  breaks = seq(from=0,to=1,length=101)
  color = colorRampPalette(c("white","red"))(100)
  pheatmap(mat = JD,scale = 'none',breaks = breaks,color = colorRampPalette(c("white","red"))(100))
}


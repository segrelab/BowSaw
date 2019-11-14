# 28 Nov 18
# Author Demetrius DiMucci

print('Loading CRAN packages')
library(randomForest)
library(pROC) 
print('Loading custom BowSaw scripts and functions')
source('synthetizeData.r')
source('BowSaw_Evaluation.r')
source('countNetwork.r')
source('findNetwork.r')
source('rankEdges.r')
source('curateNetwork.r')



### write a function that will do everything above so that we have to only specify a few parameters.
BowSawExample <- function(nGenes,nSamps,numRules,ruleLengths,rulePenetrance,overlap,
                          mixed,background,categories,trueCat,trees,minmatch,trackOOb,ask){
  if(missing(minmatch)){
    minmatch = 10
  }
  
  if(missing(trackOOb)){
    trackOOb = TRUE
  }
  
  if(missing(categories)){
    categories = 2
  }
  if(missing(trueCat)){
    trueCat = FALSE
  }
  
  ruleMat <- generateRules(nGenes = nGenes,numRules = numRules,ruleLengths = ruleLengths,overlap = overlap,mixed = mixed,categories=categories)
  
  genes <- matrix(0,ncol=nGenes,nrow=nSamps)
  # Randomly populate the matrix
  for(i in 1:nrow(genes)){
    genes[i,] <- round(runif(nGenes,min=1,max=categories))
    #genes[i,] <- rcauchy(nGenes, location = 0, scale = 1)
    # I indicate the presence with a 1 and absence with a 2 instead of 1 and 0 respectively.
    # This is just for convenience for how I've coded up BowSaw so far.
  }
  
  
  
  targetMat <- findTargets(genes = genes,ruleMat = ruleMat)
  labels <- assignLabels(targetMat = targetMat, rulePenetrance = rulePenetrance ,background = background)
  
  labs <- as.factor(labels[[1]])
  if(trueCat == TRUE){
    cats <- unique(as.vector(genes))
    for(i in 1:length(cats)){
      genes[which(genes == cats[i])] = LETTERS[i]
    }
  }
  
  print(table(labs))
  
  # if(ask == TRUE){
  #   switch(menu(choices = c('Y','N'),title = 'Is this distribution of labels OK?') ,mod <- randomForest(labs ~., genes,keep.inbag=TRUE, ntree = trees),stop('Label distribution rejected.'))
  # } else {
  #   mod <- randomForest(labs ~., genes,keep.inbag=TRUE, ntree = trees)
  # }
  mod <- randomForest(labs ~., genes,keep.inbag=TRUE, ntree = trees)
  
  # Find the initial BowSaw networks
  samps <- which(labs == 1) # Decide which labels you're interested in finding rules for
  
  BowSawNets <- matrix(0,ncol=ncol(genes),nrow=nrow(genes))
  BowSawRules <- matrix(0,ncol=ncol(genes),nrow=nrow(genes))
  buddies <- rep(0,nrow(genes))
  pb <- txtProgressBar(min = 0, max = length(samps), style = 3)
  for(i in 1:length(samps)){
    network <- countNetwork(target = samps[i],RFclass = mod,Data = genes)
    edges <- rankEdges(target = samps[i],RFclass = mod,DATA = genes,
                       Combinations = network, categorical = TRUE,threshold = 0)
    net <- findSmallNetwork(target = samps[i],edges = edges, Data = genes, 
                            RFclass = mod, real =TRUE,initRow = 1, categorical = TRUE, minmatch = minmatch)
    #buddies[i] <- net[[4]][length(net[[4]])]
    rule = sort(net[[1]])
    BowSawNets[samps[i],rule] = genes[samps[i],rule]
    BowSawRules[samps[i],rule] = 1
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  metrics <- evaluationMetrics(ruleMat = ruleMat,genes = genes,BowSawNets = BowSawNets)
  
  Rules <- unique(BowSawNets)
  ruleLens <- vector()
  for(i in 1:nrow(Rules)){
    ruleLens[i] <- length(which(Rules[i,] >0))
  }
  if(0 %in% ruleLens){
    Rules <- Rules[-which(ruleLens == 0),]
  }
  evaluated = evaluateBowSaw(RULES = Rules ,x=genes,mod = mod)
  
  results <- list()
  results[[1]] <- genes
  results[[2]] <- mod
  results[[3]] <- ruleMat
  results[[4]] <- targetMat
  results[[5]] <- labels
  results[[6]] <- BowSawRules
  results[[7]] <- metrics
  results[[8]] <- evaluated
  results[[9]] <- BowSawNets
  return(results)
  
}

# ###
# numRules <- 4
# ruleLengths <- rep(3,numRules)
# rulePenetrance <- rep(1,numRules)
# overlap = FALSE # For randomly selecting rules, do we want rules to be able to share variables?
# mixed = FALSE # If you want absence/presence combinations to matter not just presence then set this to TRUE
# 
# sim1 <- BowSawExample(nGenes = 100,nSamps = 1000,numRules = numRules,ruleLengths = ruleLengths,rulePenetrance = rulePenetrance,background = .05)
# x <- sim1[[1]]
# mod <- sim1[[2]]
# rulesMat <- sim1[[3]]
# BSN <- sim1[[6]]
# mets <- sim1[[7]]
# cand1 <- sim1[[8]][[1]]
# cand2 <- sim1[[8]][[2]]

findSubsets <- function(cand1,cand2){
  hasSubset <- rep(0,nrow(cand1))
  for(i in 1:nrow(cand1)){
    for(j in i:nrow(cand1)){
      if(j > i){
        rule1 <- which(cand1[i,]>0)
        rul1 <- cand1[i,rule1]
        rule2 <- which(cand1[j,]>0)
        rul2 <- (cand1[j,rule2])
        subset <- isSubset(rule1 = rule1,rule2 =rule2,rul1 = rul1,rul2 = rul2)
        if(subset == 1){
          #print(j)
          hasSubset[j]=1
        }
      }
      
    }
  }
  return(hasSubset)
}


#### IF YOU HAVE YOUR OWN DATA MATRIX YOU WANT TO SUBMIT WE CAN SKIP THE INITIAL SIMULATION PART
BowSaw_ExistingSample <- function(genes,mod,choice,minmatch){
  if(missing(minmatch)){
    minmatch = 1
  }
  if(missing(choice)){
    choice = 'ALL'
  }
  # Find the initial BowSaw networks
  labs = mod$y
  if(choice != 'ALL'){
    samps <- which(labs == choice) # Decide which labels you're interested in finding rules for
  } else {
    samps <- seq(nrow(genes))
  }
  
  
  BowSawNets <- matrix(0,ncol=ncol(genes),nrow=nrow(genes))
  BowSawRules <- matrix(0,ncol=ncol(genes),nrow=nrow(genes))
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  completeNet <- matrix(0,ncol=ncol(genes),nrow=ncol(genes))
  for(i in 1:length(samps)){
    network <- countNetwork(target = samps[i],RFclass = mod,Data = genes)
    completeNet <- completeNet + network
    diag(network) = 0
    edges <- rankEdges(target = samps[i],RFclass = mod,DATA = genes,
                       Combinations = network, categorical = TRUE,threshold = 0)
    net <- findSmallNetwork(target = samps[i],edges = edges, Data = genes,
                            RFclass = mod, real =TRUE,initRow = 0, categorical = TRUE,minmatch = minmatch)
    netrule = sort(net[[1]])
    BowSawNets[samps[i],netrule] = genes[samps[i],netrule]
    BowSawRules[samps[i],netrule] = 1
    setTxtProgressBar(pb, i)

  }
  close(pb)
  
  Rules <- unique(BowSawNets)
  ruleLens <- vector()
  for(i in 1:nrow(Rules)){
    ruleLens[i] <- length(which(Rules[i,] !=0))
  }
  if(0 %in% ruleLens){
    Rules <- Rules[-which(ruleLens == 0),]
  }
  evaluated = evaluateBowSaw(RULES = Rules ,x=genes,mod = mod)
  
  results <- list()
  results[[1]] <- genes
  results[[2]] <- mod
  results[[3]] <- BowSawNets
  results[[4]] <- evaluated
  results[[5]] <- BowSawRules
  results[[6]] <- completeNet
  return(results)
  
}

# ###
# numRules <- 4
# ruleLengths <- rep(3,numRules)
# rulePenetrance <- rep(1,numRules)
# overlap = FALSE # For randomly selecting rules, do we want rules to be able to share variables?
# mixed = FALSE # If you want absence/presence combinations to matter not just presence then set this to TRUE
# 
# sim1 <- BowSawExample(nGenes = 100,nSamps = 1000,numRules = numRules,ruleLengths = ruleLengths,rulePenetrance = rulePenetrance,background = .05)
# x <- sim1[[1]]
# mod <- sim1[[2]]
# rulesMat <- sim1[[3]]
# BSN <- sim1[[6]]
# mets <- sim1[[7]]
# cand1 <- sim1[[8]][[1]]
# cand2 <- sim1[[8]][[2]]
findSubsets <- function(cand1,cand2){
  hasSubset <- rep(0,nrow(cand1))
  for(i in 1:nrow(cand1)){
    for(j in i:nrow(cand1)){
      if(j > i){
        rule1 <- which(cand1[i,]>0)
        rul1 <- cand1[i,rule1]
        rule2 <- which(cand1[j,]>0)
        rul2 <- (cand1[j,rule2])
        subset <- isSubset(rule1 = rule1,rule2 =rule2,rul1 = rul1,rul2 = rul2)
        if(subset == 1){
          #print(j)
          hasSubset[j]=1
        }
      }
      
    }
  }
  return(hasSubset)
}


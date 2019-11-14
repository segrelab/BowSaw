# 28 Nov 18
# Author Demetrius DiMucci

# These functions will create a synthetic data set with labels that we can use to test
# the limits of BowSaw. The first focus will be on using binary input data - this will
# simulate presence/absence information of genes.
# This first matrix is composed of uniformly distributed elements

# If we don't have specific rules in mind we can randomly create them (this is going to be the easier method
# for when we start evaluating different scenarios)


#####################################################
#
#           BEGIN
#             DEFINING
#               GENERATING
#                 FUNCTIONS
#
#
#
#
######################################################



###############################
#
#
# generateRules()
#
#
###############################
# The function generateRules() will randomly select the variables that comprise the specified number of rules
# of the specified lengths and the values of those variables 
generateRules <- function(nGenes,numRules,overlap,mixed,ruleLengths,categories){
  if(missing(overlap)){
    overlap = FALSE
    print('overlap not specified, setting to FALSE')
  }
  if(missing(mixed)){
    mixed = FALSE
    print('mixed not specified, setting to FALSE')
  }
  
  if(missing(nGenes)){
    print("I need to know how many genes/predictors there are. Provide argument nGenes as an int.")
  }
  if(missing(numRules)){
    print("I need to know how many rules you want. Provide argment numRules as an int.")
  }
  if(missing(ruleLengths)){
    print("I need to know how long you want each rule to be. Provide argument ruleLengths as a vector of ints.")
  }
  
  ruleMat <- matrix(0,ncol=nGenes,nrow=numRules)
  if(overlap == FALSE){
    for(i in 1:nrow(ruleMat)){
      taken <- colSums(ruleMat)
      eligible <- seq(ncol(ruleMat))[which(taken == 0)]
      tempRule <- sample(eligible,ruleLengths[i])
      ruleMat[i,tempRule] = 1
      if(mixed == TRUE){
        mix <- sample(unique(categories),length(tempRule),replace = TRUE)
        ruleMat[i,tempRule] = mix
      }
    }
  } else {
    for(i in 1:nrow(ruleMat)){
      eligible <- seq(ncol(ruleMat))
      tempRule <- sample(eligible,ruleLengths[i])
      ruleMat[i,tempRule] = 1
      if(mixed == TRUE){
        mix <- sample(unique(categories),length(tempRule),replace = TRUE)
        ruleMat[i,tempRule] = mix
      }
    }
    
  }
  return(ruleMat)
}

###############################
#
#
# findTargets()
#
#
###############################
# Once we have our rules set we need to identify which samples have each rule.
# The funciton findTargets() will identify which samples in the data matrix possess the rule and
# return a matrix of targets.
findTargets <- function(genes,ruleMat){
  # Initialize a targets matrix, so we can quickly look up how many rules a given sample possesses
  # The rows will correspond to the rows of the data matrix
  # The columns will be the relevant rules - if a sample has a rule it will get marked with a '1'
  targetMat <- matrix(0,nrow=nrow(genes),ncol=nrow(ruleMat))
  
  # Identify the rules and find which samples in the data possess each rule
  for(i in 1:nrow(ruleMat)){
    rule <- which(ruleMat[i,] != 0)
    # It's not enough to know which variables are in the rule, we also need to know their values
    vals <- ruleMat[i,rule]
    
    matches <- seq(nrow(genes))
    for(j in 1:length(rule)){
      matches <- intersect(matches,which(genes[,rule[j]] == vals[j]))
    }
    targetMat[matches,i] = 1
  }
  
  return(targetMat)
}

###############################
#
#
# assignLabels()
#
#
###############################
# The last step is to assign labels based on the penetrances. Using the targetMatrix we will determine eligible
# samples and assign them a label according to the penetrance of the relevant rules.
# The function assignLabels takes the targetMatrix and rulePenetrance as arguments and will assign a label to 
# each sample from a binomial distribution equivalent to the corresponding rule penetrance.
# It will return two results - 1 the vector of labels, 2 the matrix of labels - this way one can investigate
# which rules were responsible for assigning the label if they are so inclined.
assignLabels <- function(targetMat,rulePenetrance,background){
  if(missing(background)){
    background = 0
  }
  
  labelMat <- targetMat*0
  for(i in 1:ncol(targetMat)){
    targets = which(targetMat[,i] == 1)
    prob = rulePenetrance[i]
    labels <- rbinom(length(targets),1,prob=prob)
    labelMat[targets,i] = labels
  }
  labs <- sign(rowSums(labelMat))
  
  if(background > 0 & background < 1){
    noiseLab <- rbinom(length(labs),1,prob=background)
    labs[which(noiseLab == 1)] = 1 
  }
  
  
  results <- list(); results[[1]] <- labs; results[[2]] <- labelMat
  return(results)
}

# Invoking the above 3 functions will allow us to create a data set for training a random forest and 
# subsequently applying BowSaw to.


#####################################################
#
#         END
#           DEFINING
#               GENERATING
#                   FUNCTIONS
#
#
#
#
######################################################



#####################################################
#
#           BEGIN
#             DEFINING
#               EVALUATION
#                 FUNCTIONS
#
#
#
#
######################################################


# We need some functions to evaluate the BowSaw outputs
# namely we want to determine a few things on a per-rule basis
# 1 - for each rule - how many of the correct variables were recovered for it in each network?
# 2 - for each rule - how many extraneous variables were returned in each network?
# 3 - How frequently is each variable of the rule recovered across networks?
# 4 - How many extra variables (here a variable will be extra if it doesn't belong to any rule) are found across all networks?

# Find how many correct variable+values from each rule were found in each BowSaw network.
# And records how many noise elements were returned
evaluationMetrics <- function(ruleMat,genes,BowSawNets){
  foundMat <- matrix(0,nrow=nrow(genes),ncol=nrow(ruleMat))
  haveMat <- foundMat
  extraMat <- foundMat
  foundVectors <- ruleMat
  extraVectors <- matrix(0,ncol(genes),nrow=length(unique(as.vector(genes))))
  
  allRules <- which(ruleMat != 0, arr.ind = TRUE)
  allVals <- ruleMat[allRules]
  allVars <- allRules[,2]
  
  
  for(i in 1:nrow(ruleMat)){
    rule <- which(ruleMat[i,] != 0)
    vals <- ruleMat[i,rule]
    
    for(j in 1:nrow(BowSawNets)){
      foundMat[j,i] = length(which(BowSawNets[j,rule] == vals)) # answers metric 1
      extraMat[j,i] = length(which(BowSawNets[j,-rule] != 0)) # answers metric 2
      haveMat[j,i] = length(which(genes[j,rule] == vals))
    }
    
    # answers metric 3
    found <- rep(0,length(rule))
    for(j in 1:length(rule)){
      found[j] <- length(which(BowSawNets[,rule[j]] == vals[j]))
    }
    foundVectors[i,rule] <- found
    
    # answers metric 4 -  records which extra variables+values were found and how frequently
    # variable 4 of category 2 is not the same as variable 4 of category 1.
    # we should therefore differentiate between categories
    categories <- sort(unique(as.vector(genes)))
    for(j in 1:ncol(BowSawNets)){
      if(!(j %in% allVars)){
        for(k in 1:length(categories)){
          extraVectors[k,j] <- length(which(BowSawNets[,j] == categories[k]))
        }
        
      } else {
        relevant <- which(allVars == j)
        for(k in length(categories)){
          extraVectors[k,j] <- length(which(BowSawNets[,j] == categories[k] & !(BowSawNets[,j] %in% allVals[relevant])))
        }
        
      }
    }
  }
  
  results <- list()
  results[[1]] <- foundMat # for each BowSaw network, for each rule how many variables+values of that rule were returned
  results[[2]] <- extraMat # for each BowSaw network, for each rule how many variables+values that were NOT that rule were returned
  results[[3]] <- foundVectors # for each rule how many times each variable+value was returned
  results[[4]] <- extraVectors # how many times a variable+value that is not in a rule was returned
  results[[5]] <- allVars # Which variables are parts of rules
  results[[6]] <- allVals # what the corresponding values/categories are
  results[[7]] <- haveMat # For each sample, for each rule how many variables+values of that rule are present
  return(results)
}

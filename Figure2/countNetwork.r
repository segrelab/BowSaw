# July 26 
# Author: Demetrius DiMucci

# This function will, for a single observation: 
# Find all proximate pairwise combinations (i -> j) in out of bag trees

# 3 arguments must be provided:
# 1 - target: this is the row index of the sample you're interested in
# 2 - RFclass: this is the random forest model you have trained
# 3 - Data: this is the feature matrix used to train the model

countNetwork <- function(target,RFclass,Data){
  
  
  Forest <- RFclass$forest
  OOb <- RFclass$oob.times
  inbag <- RFclass$inbag
  Labels <- RFclass$y
  predictLabels <- RFclass$predicted
  
  classes <- unique(RFclass$y)
  
  # Identify for which trees in the forest the sample of interest was out of bag.
  # These are the trees we will evaluate
  # You can modify this to look at the inbag trees if you wish. Also the argument trackOOb, if provided will switch you to
  # looking at ALL tree paths
  oobSamps <- which(inbag[target,]==0) 

  # Extract the feature vector of the sample
  # We use this to find the paths of our sample
  Case <- Data[target,]

  # Keep track of how often each combination occurs. Note that a V6 -> V3 split is different from V3 -> V6 split.
  combinations <- matrix(0,ncol = ncol(Data),nrow=ncol(Data))

  # Determine what the mean probability of being labeled class '1' is
  # if we just used the distribution of the full data set to guess.
  #####
  for(i in 1:length(oobSamps)){
    # Extract an individual tree from the forest
    Tree <- getTree(RFclass,oobSamps[i],labelVar = F)
    status = 1 # Status - means it's a terminal node
    node = 1 # We start at the root node, i.e. node 1
    path <- 0 # No decision has been made yet while in the root node
    vars <- 0 # ditto
    while(status == 1){
      left <- 0 # use this variable to determine which node is next in the path
      splitVar <- Tree[node,3] # Which variable was used to split
      
      path <-c(path,node) # Keep track of the nodes the sample visited
      
      vars <- c(vars,splitVar) # Keep track of the decisions made at each node
      
      splitPoint <- Tree[node,4] # identify the decision variable at the current node
      
      if(Case[splitVar] <= splitPoint){
        left = 1
        node <- Tree[node,1]
      } else{
        left = 0
        node <- Tree[node,2]
      }
      # Update the node type, if it is terminal we will stop
      status <- Tree[node,5] # -1 is a terminal node
      
      # Dynamically update the combinations matrix
      if(length(path) > 1){
        h <- length(path)
        # Count how many times each combination shows up for this sample
        # Columns are the previous splitting variable
        # Rows are the current splitting variable
        combinations[vars[h],vars[h-1]]=combinations[vars[h],vars[h-1]]+1
      }
    }
  }

  results <- combinations #These are immediate j -> i edges.
  return(results)
}

# Test it 
# target=3
# S <- countNetwork(target,RFclass,Data)
# network <- countNetwork(target = samps[i],RFclass = mod,Data = genes)

# 7 March 2019
# Demetrius DiMucci
# 


# FUNCTION THAT FINDS ALL POSSIBLE KMERS THAT CAN BE BUILT FROM BOWSAW RULES
kmerFind <- function(candList,kmers){
  kmerlist <- apply(candList,1,function(x) combn(paste0(which(x>0),"_",x[which(x>0)]),kmers))
  kmer1 <- do.call(cbind,kmerlist)
  
  k1=apply(format(kmer1), 2, paste0, collapse="_")
  k2 = gsub(pattern = ' ',replacement = '',x = k1)
  k3 = table(k2)
  
  return(k3)
}

## EVLAUTE THE KMERS
evaluateKmer <- function(t11,x,labs,rulesMat){
  fullSizes <- rowSums(sign(rulesMat))
  matchVec <- vector()
  LabVec <- vector()
  t11Index <- vector()
  
  valSet <- vector()
  ruleSet <- vector()
  
  fullMatch <- rep(0,length(t11))
  ruleNames <- vector()
  ruleVals <- vector()
  
  exhaustCand <- matrix(0,nrow=length(t11),ncol=ncol(x))
  ################
  
  
  for(i in 1:length(t11)){
    # BREAK UP THE RULE FROM T11 INTO THE RULE AND THE CORRESPONDING VALUES
    rule_val <- as.numeric(strsplit(names(t11)[i],"_")[[1]])
    rule <- rule_val[seq(1,length(rule_val),2)]
    vals <- rule_val[seq(1,length(rule_val),2)+1]
    
    exhaustCand[i,rule] = vals
    
    ruleNames[i] <- paste(rule,collapse = '_')
    ruleVals[i] <- paste(vals,collapse = '_')
    
    # HOW MANY VARIABLES OF EACH REAL RULE ARE IN THIS RULE?
    sizes <- rowSums(sign(rulesMat[,rule]))
    
    # DETERMINE IF ANY OF THE REAL RULES LOOK LIKE THEY'RE CONTAINED IN THIS RULE
    hit1 <- fullSizes - sizes
    
    if(0 %in% hit1){
      # IF TRUE THEN THE VARIABLES OF AT LEAST ONE REAL RULE ARE IN THIS CANDIDATE RULE
      hit = which(hit1 == 0)
      
      # NEXT CHECK TO SEE IF THE VALUES OF THE VARIABLES MATCH UP TO THE POTENTIAL MATCH
      tempRule <- which(rulesMat[hit,]>0)
      lilRule = rule[which(rule %in% tempRule)]
      lilVal = vals[which(rule %in% tempRule)]
      tempVal = rulesMat[hit,lilRule]
      
      if(sum(abs(tempVal - lilVal))==0){
        # IF THE VARIABLES MATCH UP TO ONE OF THE RULES RECORD WHICH RULE IS CONTAINED IN THE CANDIDATE RULE
        fullMatch[i]=hit
        if(length(tempRule) == 4){
          fullMatch[i] = hit 
        }
        # IF ITS AN EXACT MATCH MARK IT 
        if(length(rule) == length(tempRule)){
          fullMatch[i] = hit + nrow(rulesMat)
        }
        
      } else {
        fullMatch[i] = 0
      }
    }
    
    # NEXT, RECORD HOW MANY SAMPLES IN THE DATA HAVE THIS RULE
    x2 <- x[,rule]
    z <- apply(x2, 1, function(a) all(a == vals))
    keep <- which(z == TRUE)
    
    matchVec[i] <- length(keep)
    # RECORD HOW MANY SAMPLES HAVE THE RULE AND HAVE THE LABEL
    LabVec[i] <- length(which(labs[keep] == 1))
  }
  
  
  results <- list()
  results[[1]] <- fullMatch
  results[[2]] <- matchVec
  results[[3]] <- LabVec
  results[[4]] <- exhaustCand
  
  return(results)
}





# ### test it out for something
# t11 = kmerFind(candList = cand1,kmers = 4)
# 
# t1 <- evaluateKmer(t11 = t11,x = x, labs = labs,rulesMat = rulesMat)
# 
# 
# mers <- seq(5)[-1]
# 
# finals <- list()
# for(i in 1:length(mers)){
#   t11 <- kmerFind(cand3,kmers=mers[i])
#   finals[[i]] <- evaluateKmer(t11 = t11, x=x,labs = labs, rulesMat = rulesMat)
# }
# 
# plot(finals[[4]][[2]],finals[[4]][[3]]/finals[[4]][[2]])
# 
# 
# 
# 
# allMatch <- vector()
# allLabeled <- vector()
# allCol <- vector()
# for(i in 1:length(mers)){
#   allMatch <- c(allMatch,finals[[i]][[2]])
#   allLabeled <- c(allLabeled,finals[[i]][[3]])
#   allCol <- c(allCol,finals[[i]][[1]])
# }
# 
# colors = c('white','gray','black','blue','red','darkgreen','brown')
# colors2 = allCol + 1
# 
# plot(allMatch,allLabeled/allMatch,xlab='Matching Samples',ylab='Fraction With Label')
# 
# d=allLabeled/allMatch
# 
# keep <- which(colors2 != 1)
# points(allMatch[keep],d[keep],col=colors2[keep],pch=16)

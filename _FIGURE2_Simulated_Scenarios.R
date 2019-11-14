# 8 March 2019
# Run A simulation with parameters that give a similar output to what We seee with the drosophila enhancer data.

library(pROC)
library(randomForest)
library(pheatmap)
library(PRROC)
source('Functions.r')
source('prioritize.r')
source('kMerFinding.r')

set.seed(5)

numRules <- 5 # HOW MANY WAYS THERE ARE TO GET THE LABEL/PHENOTYPE
# Simulation 1
ruleLengths <- round(runif(numRules,min=4,max=4)) # HOW COMPLICATED THE RULES ARE
NOISE = 0 # THE PROBABILTY OF SPONTANEOUSLY HAVING THE LABEL REGARDLESS OF RULE STATUS

# # Simulation 2
# Uncomment this to get the more challenging scenario
# ruleLengths <- round(runif(numRules,min=3,max=5)) # HOW COMPLICATED THE RULES ARE
# NOISE = 0.02 # THE PROBABILTY OF SPONTANEOUSLY HAVING THE LABEL REGARDLESS OF RULE STATUS

rulePenetrance <- runif(numRules,min=.8,max=.9)  # PROBABILTY OF GETTING THE LABEL GIVEN THE RULE
overlap = TRUE # For randomly selecting rules, do we want rules to be able to share variables?
mixed = TRUE # IF SET TO FALSE ALL VARIABLES VALUES MUST BE 1.
#IF SET TO TRUE VARIABLE VALUES WILL BE RANDOMLY SELECTED FROM THE NUMBER OF CATEGORIES

CATS = 2 # HOW MANY VARIABLE CATEGORIES ARE THERE?
totTrees = 500 # HOW MANY TREES YOU WANT TO CREATE
samples = 2000




sim2 <- BowSawExample(nGenes = 100,nSamps = samples,numRules = numRules,trees = totTrees,overlap = overlap,mixed=mixed,
                      ruleLengths = ruleLengths,rulePenetrance = rulePenetrance,background = NOISE,
                      categories = CATS,minmatch = 1,trackOOb = TRUE,ask=F)


x <- sim2[[1]]
mod <- sim2[[2]]
labs <- mod$y
rulesMat <- sim2[[3]]
roc(mod$y,mod$votes[,2])
cand1 <- sim2[[8]][[1]] # The matrix of unique BowSaw networks identified
cand2 <- sim2[[8]][[2]] # Initial evaluation of each BowSaw network sorted in descending order of how many observed samples the rule captures


varImpPlot(mod,main='')

rez = prioritize(cand1 = cand1,cand2 = cand2,x = x,mod = mod,lab = 1)
rezO <- unchanged(cand1 = cand1,cand2=cand2,x=x,mod=mod,lab=1)
accounted <- rezO[[4]]

newCand1a <- rez[[1]]
newCand2a <- rez[[2]]
accounted1 <- rez[[4]]
bor <- rez[[3]]
added <- rez[[5]]
last <- min(which(accounted1== 0))-1
newCand1a = cand1
labs = mod$y

plot(accounted,type='l',lwd=3,ylab='Fraction Accounted For',xlab='Candidate Rule Rank')
points(accounted1[1:last],type='l',lwd=3,col=2)
legend('bottomright',legend=c('Original','Reprioritized'),pch=16,col=c(1,2))

cand3 <- cand1[bor,]
cand4 <- cand2[bor,]

########## FIND THE INDICES OF THE FIRST FULL APPEARANCE OF EACH OF THE 5 RULES

index1 <- rep(0,numRules)
index2 <- index1

for(i in 1:nrow(rulesMat)){
  rule = which(rulesMat[i,]>0)
  vals = rulesMat[i,rule]
  
  for(j in 1:nrow(cand1)){
    g = all.equal(cand1[j,rule],vals)
    if(g == TRUE & index1[i] == 0){
      index1[i] = j
    }
  }
  
  for(j in 1:nrow(cand3)){
    g = all.equal(cand3[j,rule],vals)
    if(g == TRUE & index2[i] == 0){
      index2[i] = j
    }
  }
}
index1
index2


########### 
ts <- vector()
potential <- 0
for(i in 2:5){
  
  small <- which(cand2[,4] < i)
  if(length(small) > 0){
    tempCand = cand3[-small,]
  } else {
    tempCand = cand3
  }
  
  ts <- c(ts,kmerFind(candList = tempCand,kmers = i))
  
  
  potential <- potential + choose(100,i)*2^i
}
length(ts)

# APPLY ALGORITHM 4
p=evaluateKmer(t11 = ts, x=x,labs = labs, rulesMat = rulesMat)

b = p[[2]]
c = p[[3]]
d = c/b

keep <- which(p[[1]]==0)
keep1 <- which(p[[1]] > 0)
keep2 <- which(p[[1]] > 0 )

pdf()
plot(b[keep],d[keep],ylim=c(0,1),xlab='Matching Samples',ylab='Precision',cex.lab=1.5,cex.axis=1.2)
colahs = c('red','black','brown','brown','darkgreen')
points(b[keep2],d[keep2],col=1)
keep3 <- which(p[[1]] > 5)
points(b[keep3],d[keep3],col=2,pch=16)
thresh = .75
perfs <- which(d >= thresh)
abline(h=thresh,lty=3)
legend('topright',legend=c('Candidate Rule','True Rule'),pch=16,col=c(1,2),cex=1.2)
dev.off()

######


###### APPLY THE MODIFIED VERSION OF PRIORITIZATION ALGORITHM 3
#### ALGORITHM 5!



bestboi = perfs[which.max(b[perfs])]
roolz = c(ts[bestboi],ts)

tCand <- p[[4]]
tCand2 <- cbind(b,d,c)
tCand <- rbind(tCand[bestboi,],tCand)
tCand2 <- rbind(tCand2[bestboi,],tCand2)
consider <- which(tCand2[,2] >= thresh)
roolz3 <- roolz[consider]
tCand1 <- tCand[consider,]
tCand3 <- tCand2[consider,]

tRez <- prioritize_imperfect(cand1 = tCand1,cand2=tCand3,x=x,mod=mod,lab=1)


accounted1 <- tRez[[4]]
bor <- tRez[[3]]
added <- tRez[[5]]

keep <- max(which(bor != 0))
discovered <- rep(0,(keep))
sizes <- rowSums(sign(rulesMat))
for(i in 1:keep){
  rule_val <- as.numeric(strsplit(names(roolz3)[bor[i]],"_")[[1]])
  rule <- rule_val[seq(1,length(rule_val),2)]
  vals <- rule_val[seq(1,length(rule_val),2)+1]
  
  ### check if the rule is in the rulesMat - if it is an exact match give a 1 if not give a 0
  for(j in 1:nrow(rulesMat)){
    tempRule <- which(rulesMat[j,]>0)
    g = all.equal(tempRule,rule)
    if(g == TRUE){
      h = all.equal(rulesMat[j,tempRule],vals)
      if(h == TRUE){
        discovered[i]=j
      }
    }
  }
}
discovered



precision <- cumsum(sign(discovered))/seq(length(discovered))
recall <- cumsum(sign(discovered))/5
plot(recall,precision,type='l')

precision <- c(1,precision,1)
recall <- c(0,recall,1)
plot(recall,precision,type='l')
auc(recall,precision)

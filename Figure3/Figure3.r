# 9 Feb 2019
# Script for making the panels for Figure 5
# Author: Demetrius DiMucci

# THE DATA FOR THIS PART CAME FROM HERE: https://ibdmdb.org/tunnel/public/summary.html
setwd('C:/Users/ddimucci/Desktop/Manuscripts/BowSaw/Figure3/')

library(pROC)
library(randomForest)
#library(pheatmap)
library(PRROC)
#library(precrec)
source('Functions.r')
source('prioritize.r')
source('kMerFinding.r')


##### LOAD THE DATA AND DO A BIT OF CLEANING TO BE ABLE TO USE IT
 x <- read.csv('hmp2_metadata.csv')
 y <- read.csv('taxonomic_profiles.tsv',sep='\t')

phyl = y[,ncol(y)]
OTU <- y[,1]

row.names(y) = y[,1]
y = y[,-1]
y = y[,-ncol(y)]

IDs <- names(y)
IDs <- gsub('X','',IDs)
colnames(y) = IDs

OTUtable <- t(y)
PRESENCEtable <- sign(t(y))
# MAKE THE VALUES THAT ARE SET TO ZERO A NON-ZERO NUMBER BECAUSE THE CODE FOR BOWSAW USES 0 TO INDICATE A VARIABLE IS TO BE IGNORED
PRESENCEtable[which(PRESENCEtable == 0)] = 2 
# THERE ARE MORE DIAGNOSES THAN OTU PROFILES - IDENTIFY WHICH DIAGNOSES SAMPLES ALSO HAVE OTU PROFILES AND KEEP THEM
ord <- vector()
for(i in 1:length(IDs)){
  ord[i] <- which(x$External.ID == IDs[i])
}

diagnosis <- x$diagnosis[ord]
diagnosis <- as.factor(diagnosis)

location <- x$biopsy_location[ord]


cdonly <- rep(0,length(diagnosis))
cdonly[which(diagnosis == 'CD')]=1
cdonly <- as.factor(cdonly)

CDmod <- randomForest(cdonly ~., OTUtable,keep.inbag=TRUE,ntree=2000)
CDmodBinary <- randomForest(cdonly ~., PRESENCEtable,keep.inbag = TRUE,ntree=2000)

pred <- CDmod$y
votes = CDmod$votes[,2]
g1 = which(pred == 1)
fg1 = votes[g1]
bg1 = votes[-g1]
#auc1 = evalmod(scores = CDmod$votes[,2],labels=CDmod$y)
roc1 = roc.curve(scores.class0 = fg1,scores.class1 = bg1,curve = TRUE)
pr1 = pr.curve(scores.class0 = fg1,scores.class1 = bg1,curve = TRUE)


pred <- CDmodBinary$y
votes = CDmodBinary$votes[,2]
fg1 = votes[g1]
bg1 = votes[-g1]
roc2 = roc.curve(scores.class0 = fg1,scores.class1 = bg1,curve = TRUE)
pr2 = pr.curve(scores.class0 = fg1,scores.class1 = bg1,curve = TRUE)

# Panel A
jpeg(file='ROCs.jpeg')
plot(roc1,legend = F,col='black',auc.main = F,main = paste('ROC Curves'),cex.lab=1.5,cex.main=1.5,xlab='1 - Specificity')
plot(roc2,add = T,col='red')

aucs <- vector(); aucs[1] <- paste('Original:',round(roc1$auc,3));aucs[2] <- paste('Discretized:',round(roc2$auc,3))
legend('bottomright',legend=aucs,pch=16,col=c(1,2),title='AUC')
abline(0,1,lty=3)
dev.off()

# Panel B
jpeg(file='PRS.jpeg')
plot(pr1,legend = F,col='black',auc.main = F,main = paste('PR Curves'),cex.lab=1.5,cex.main=1.5)
plot(pr2,add = T,col='red')

aucs <- vector(); aucs[1] <- paste('Original:',round(pr1$auc.integral,3));aucs[2] <- paste('Discretized:',round(pr2$auc.integral,3))
abline(h=86/nrow(PRESENCEtable),lty=3)
legend('bottomright',legend=aucs,pch=16,col=c(1,2),title='AUC')
dev.off()

### MAKE PANEL C
### HOW MANY RULES DO WE NEED TO ACCOUNT FOR EVERYTHING?

# cdBow <- BowSaw_ExistingSample(genes = PRESENCEtable,mod = CDmodBinary,choice = 1,minmatch = 1)
# 
# save(cdBow,file='cdBow')
# Load a pre-computed BowSaw results
load('cdBow')
PRESENCEtable <- cdBow[[1]]
CDmodBinary <- cdBow[[2]]
cand1 <- cdBow[[4]][[1]]
cand2 <- cdBow[[4]][[2]]


rez <- prioritize(cand1 = cand1,cand2=cand2,x=PRESENCEtable,mod=CDmodBinary,lab=1)
rezO <- unchanged(cand1 = cand1,cand2=cand2,x=PRESENCEtable,mod=CDmodBinary,lab=1)

newCand1 <- rez[[1]]
newCand2 <- rez[[2]]
bor = rez[[3]]
accounted <- rez[[4]]
added <- rez[[5]]

last <- min(which(accounted== 0)) - 1
plot(accounted[1:last],xlab='Candidate Rule Index',ylim=c(0,1),ylab='Fraction of Samples Accounted For')

# Exhaustively examine the sub-rules
cand3 <- newCand1[1:last,]
cand2 <- cand2[bor,]

bugRules <- vector()
bugLengths = vector()
for(i in 2:7){
  
  small <- which(cand2[,4] < i)
  
  if(length(small) > 0){
    tempCand = cand3[-small,]
  } else {
    tempCand = cand3
  }
  
  bugRules <- c(bugRules,kmerFind(candList = tempCand,kmers = i))
  j = kmerFind(candList = tempCand,kmers = i)
  bugLengths = c(bugLengths,rep(i,length(j)))
}
length(bugRules)



####################
fakeMat <- matrix(1,nrow=5,ncol=982) # Just feed a fake rules matrix to the function when using real data
cdonly <- CDmodBinary$y
CDrules =evaluateKmer(t11 = bugRules, x=PRESENCEtable,labs = cdonly, rulesMat = fakeMat)

eCand <- CDrules[[4]]

b1 = CDrules[[2]]
c2 = CDrules[[3]]
d2 = c2/b1

# PANEL C
plot(b1,d2,xlab='Number of Exactly Matching Samples',ylab="Precision",cex.lab=1.5)
abline(h=.9)

perfs <- which(d2 >= 1)
bestboi = perfs[which.max(b1[perfs])]
bug2 = c(bugRules[bestboi],bugRules)

eCand2 <- cbind(b1,d2,c2)
eCand <- rbind(eCand[bestboi,],eCand)
eCand2 <- rbind(eCand2[bestboi,],eCand2)



consider <- which(eCand2[,2] >= .9)
bug3 <- bug2[consider]
eCand1 <- eCand[consider,]
eCand3 <- eCand2[consider,]
rez <- prioritize_imperfect(cand1 = eCand1,cand2=eCand3,x=PRESENCEtable,mod=CDmodBinary,lab=1)

best <- rez[[3]]

accounted <- rez[[4]]
bor <- rez[[3]]
added <- rez[[5]]
bug3[best]
#accounted
t2 = bug3[best]
j=1
rule_val <- as.numeric(strsplit(names(t2)[j],"_")[[1]])
rule <- rule_val[seq(1,length(rule_val),2)]
vals <- rule_val[seq(1,length(rule_val),2)+1]

leMatch <- seq(nrow(PRESENCEtable))
for(i in 1:length(rule)){
  leMatch <- intersect(leMatch,which(PRESENCEtable[,rule[i]]==vals[i]))
}
table(cdonly[leMatch])

phyl[rule]



### Display the rules
for(i in 1:length(best)){
  j=i
  rule_val <- as.numeric(strsplit(names(t2)[j],"_")[[1]])
  rule <- rule_val[seq(1,length(rule_val),2)]
  vals <- rule_val[seq(1,length(rule_val),2)+1]
  
  leMatch <- seq(nrow(PRESENCEtable))
  for(i in 1:length(rule)){
    leMatch <- intersect(leMatch,which(PRESENCEtable[,rule[i]]==vals[i]))
  }
  tab = table(cdonly[leMatch])
  
  phyl[rule]
  print(phyl[rule])
  print(c('y','n')[vals])
  print(paste('************',j,tab[1],tab[2]))
}

########### Check how the aggregate rules behave
# Obtain results for table 2
# I Manually wrote these into the table
jpeg()
plot(b1,d2,xlab='Number of Exactly Matching Samples',ylab="Precision",cex.lab=1.5)
abline(h=.9)
allMatch <- vector()
for(j in 1:length(t2)){
  rule_val <- as.numeric(strsplit(names(t2)[j],"_")[[1]])
  rule <- rule_val[seq(1,length(rule_val),2)]
  vals <- rule_val[seq(1,length(rule_val),2)+1]
  
  leMatch <- seq(nrow(PRESENCEtable))
  for(i in 1:length(rule)){
    leMatch <- intersect(leMatch,which(PRESENCEtable[,rule[i]]==vals[i]))
  }
  
  hits = length(leMatch)
  prec = length(which(CDmodBinary$y[leMatch] == 1))/hits
  
  points(hits,prec,col=2,pch=16,cex=1.5)
  legend('bottomright',legend=c('Selected Sub-Rule'),col=2,pch=16)
  print(table(cdonly[leMatch]))
  #print(phyl[rule])
  #print(table(cdonly[allMatch]))
  allMatch <- union(allMatch,leMatch)
}
dev.off()
table(cdonly[allMatch])
table(diagnosis[allMatch])

labs2 = rep(0,length(cdonly))
labs2[allMatch]=1
auc(cdonly,labs2)


##########
isparent = vector()
parents = rowSums(sign(cand3))
parlength = vector()
chillength = vector()
for(i in 1:length(best)){
  r = names(bug3[best[i]])
  r = gsub(pattern = '_1_',replacement = ' ',x = r)
  r = gsub(pattern = '_2_',replacement = ' ',x = r)
  r = gsub(pattern = '_1',replacement = ' ',x = r)
  r = gsub(pattern = '_2',replacement = ' ',x = r)
  r = as.numeric(strsplit(x = r,split = ' ')[[1]])
  
  childof = which(rowSums(sign(cand3[,r])) == length(r))
  print(childof)
  parlength[i] = parents[childof]
  chillength[i] = length(r)
  isparent[i] = length(r) %in% parents[childof]
}
summary(chillength)
summary(parlength)




################## If you're curious about how logistic regression does on this data set
B = (OTUtable)
labs = CDmod$y


train = sample(nrow(B),90)
B1 = data.frame(B[train,])
L1 = labs[train]
B2 = B[-train,]
train.df=data.frame(y=L1,B1)
mood=glm(y~.,train.df,family=binomial('logit'))
test.df = data.frame(y=labs[-train],B2)
pred1 = predict(mood,test.df)
plot(evalmod(scores=pred1,labels = labs[-train]))
table(sign(pred1),labs[-train])

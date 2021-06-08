# 10 June 2019
# Script for making the panels for Figure S2
# Author: Demetrius DiMucci
setwd('C:/Users/ddimucci/Desktop/Manuscripts/BowSaw/Figure_S1/')
library(randomForest)
library(plyr)


##### LOAD THE DATA FOR THE SIMULATIONS WITH 2000 SAMPLES, 200 FEATURES, 4-5 RULES 
######## VARIABLE PENETRANCE AND NOISE FROM 0 TO .45
x <- read.csv('netResults_fromCommandline.csv')
x <- x[,-1]

##### PANEL A - Show that binarizing the data results in a model that is effectively just as 
# Good as using the continuous data
columnNames <- c('minPen','maxPen','backNoise','ROC-AUC','PR-AUC','nBowSawPairs','BowSawPairPR-AUC','ExhaustivePairPR-AUC',
                 'nCandRules','Pen1','Pen2','Pen3','Pen4','Length1','Length2','Length3','Length4','Found1','Found2','Found3','Found4',
                 'Index1','Index2','Index3','Index4','NewIndex1','NewIndex2','NewIndex3','NewIndex4')
names(x) <- columnNames

### Average length of the rules
avgLength <- rowMeans(x[,14:17])

#### DETERMINE THE AVERAGE BEST RECOVERY OF THE FOUR RULES 

fractions <- x[,18:21]/x[,14:17]
meanFrac <- rowMeans(fractions)

fracs <- fractions
fracs[which(fracs < 1,arr.ind = T)]=0
completes1 <- rowSums(fracs)
######### DETERMINE THE FURTHEST DOWN THE CANDIDATE RULES NEEDED TO INVESTIGATE BEFORE ENCOUNTERING THE LAST OF THE BEST RECOVERY
maxDepth <- vector()
for(i in 1:nrow(x)){
  maxDepth[i] <- max(x[i,26:29])
}

noises = sort(unique(x$backNoise))
############## EXAMINE THE EFFECT OF BACKGROUND NOISE ON MODEL PERFORMANCE
plot(x$backNoise,x$`ROC-AUC`,xlab='Background Noise',ylab='ROC-AUC')


noiseRoc <- vector()
for(i in 1:length(noises)){
  noiseRoc[i] <- mean(x$`ROC-AUC`[which(x$backNoise==noises[i])])
}
lines(noises,noiseRoc,lwd=3,col='red')

######## PANEL B
############## EXAMINE THE RELATIONSHIP OF ROC-AUC AND NOISE TO BOWSAW PAIR PR-AUC
plot(x$backNoise,x$`BowSawPairPR-AUC`,xlab='Background Noise',ylab='Pairwise PR-AUC')

noisePR <- vector()
for(i in 1:length(noises)){
  noisePR[i] <- mean(x$`BowSawPairPR-AUC`[which(x$backNoise == noises[i])])
}
lines(noises,noisePR,lwd=3,col='red')

plot(x$`ROC-AUC`,x$`BowSawPairPR-AUC`,xlab='ROC-AUC',ylab = 'Pairwise PR-AUC')

##### PANEL C
############## EXAMINE THE RELATIONSHIP OF ROC-AUC AND NOISE TO THE AVERAGE FRACTION RECOVERED
plot(x$backNoise,meanFrac,xlab='Background Noise',ylab='Average Fraction Recovered')
noiseFrac <- vector()
for(i in 1:length(noises)){
  noiseFrac[i] <- mean(meanFrac[which(x$backNoise == noises[i])])
}
lines(noises,noiseFrac,lwd=3,col='red')

plot(x$`ROC-AUC`,meanFrac,xlab='ROC-AUC',ylab='Average Fraction Recovered')
plot(x$`ROC-AUC`,maxDepth,xlab='ROC-AUC',ylab='Candidate Rules Evaluated')

###################
########## PANEL D/E - IF WE MANAGE TO HAVE A RELATIVELY LOW-NOISE SAMPLING HOW DO THE OTHER UNDERLYING PARAMETS AFFECT 
########## MODEL QUALITY AND BOWSAW RECOVERY?
values <- c('minPen','maxPen','noise','feats','numRules','len1','len2','len3','len4','len5','len6','len7','len8',
            'pen1','pen2','pen3','pen4','pen5','pen6','pen7','pen8','coin1','coin2','ROCAUC','PRAUC','BSpairs','PosPairs','PAIR-PRAUC')
z <- read.csv('MassiveResults.csv')[,-1]
y <- read.csv('MassiveResults_small.csv')[,-1]

z <- rbind(z,y)
size <- rep(0,20000);size[1:10000]=1 # THERE WERE 2000 SAMPLES IN EACH SIMULATION OF THE Z MATRIX AND 200 IN THE Y MATRIX
z <- cbind(z,size)
colnames(z)[1:28] <- values

means <- vector() # average fraction of a rule discovered
deepest <- vector()
avgrank = vector()
rankCompletes = vector()
bestrank = vector()
completes <- vector() # how many complete rules discovered
fracCompletes <- vector() # fraction of all rules that are complete
lengths <- vector()
pens <- vector()
for(i in 1:nrow(z)){
  ord <- which(z[i,45:52] > -1) # how many variables of each rule were recovered by the best bowsaw rule
  rules <- which(z[i,6:13] > -1) # How many rules were in the simulation
  means[i] <- mean(as.numeric(z[i,45:52][ord]/z[i,6:13][ord])) # average fraction of the true rules were recovered 
  completes[i] <- length(which(as.numeric(z[i,45:52][ord]/z[i,6:13][ord]) == 1)) # how many rules were fully recovered at least once
  completed = (which(as.numeric(z[i,45:52][ord]/z[i,6:13][ord]) == 1))
  ord2 = ord[completed]
  fracCompletes[i] <- completes[i]/length(rules) 
  deepest[i] <- max(as.numeric(z[i,53:60][ord])) # by ranking the extracted bowsaw rule what is the worst rank of any of these best recovery rules
  avgrank[i] <- mean(as.numeric(z[i,53:60][ord])) # what is the average rank of these recovery rules (relative to all rules extracted)
  lengths[i] <- mean(as.numeric(z[i,6:13][ord]))
  pens[i] <- mean(as.numeric(z[i,14:21][ord]))
  if(length(ord2) > 0){
    rankCompletes[i] = median(as.numeric(z[i,53:60][ord2]))
    bestrank[i] = min(as.numeric(z[i,53:60][ord2]))
  } else {
    bestrank[i] = 20000
    if(i <= nrow(y)){
      rankCompletes[i] = 20000
    } else {
      rankCompletes[i] = 20000
    }
  }
}
nofind = which(rankCompletes > 2000)
summary(z$ROCAUC[nofind])
summary(z$ROCAUC[-nofind])

complete2 <- c(completes,completes1)
allLengths <- c(lengths,avgLength)
theDeep = c(deepest,maxDepth)

allmeans <- c(means,meanFrac)
allPRB <- c(z[,28],x$`BowSawPairPR-AUC`)
allROC <- c(z[,24],x$`ROC-AUC`)
allPR <- c(z[,25],x$`PR-AUC`)
tots <- nrow(z)+nrow(x)
M <- matrix(0,nrow=tots,ncol=11)
for(i in 1:3){
  M[,i] <- c(z[,i],x[,i])
}

M[,4] <- c(z[,5],rep(4,nrow(x)))
M[,5] <- c(z[,4],rep(200,nrow(x)))
M[,6] <- c(z[,24],x[,4])
M[,7] <- c(z[,25],x[,5])
M[,8] <- c(z[,26],x[,6])
M[,9] <- c(z[,69],rep(1,nrow(x)))
M[,10] <- c(z[,28],x$`BowSawPairPR-AUC`)
M[,11] <- c(z[,27],rep(79600,nrow(x)))

names(M) = c('min','max','noise','nRules','features','ROCAUC','PRAUC','BowPairs','Observations','BowPairAUC')
M <- as.data.frame(M)

############
train = seq(nrow(M))
train = seq(20000)
mins <- M[train,1]
maxes <- M[train,2]
noises <- M[train,3]
nrules <- M[train,4]
features <- M[train,5]
ROCAUC <- M[train,6]
PRAUC <- M[train,7]
Bpairs <- M[train,8]
sizes <- M[train,9]
posPairs <- M[train,11]



### Relationship between the average fraction of a rule recovered per candidate rule and the observable 
### Elements of the data set (number of features, ROC-AUC, PR-AUC, and number of samples)
inputs <- z[,c(4,24,25)]
inputs <- as.matrix(inputs); inputs <- cbind(inputs,sizes)
lm1 <- lm(allmeans[train] ~ inputs)
corr = round(cor(lm1$fitted.values,allmeans[train]),3)
plot(allmeans[train],lm1$fitted.values,xlab='Mean Fraction of Rules Recovered',ylab='Fitted Fraction of Rules Recovered',ylim=range(allmeans))

abline(0,1,lwd=3,col=2)
#legend('topleft',#legend=paste('r: ',corr))

# 
inputs <- z[,c(1:5)]

for(i in 1:ncol(inputs)){
  bad <- which(inputs[,i]==-1)
  inputs[bad,i]=0
}

good <- vector()
pens <- vector()
lengths <- vector()
for(i in 1:nrow(z)){
  good <- which(z[i,6:13] > 0)
  lengths[i] <- mean(as.numeric(z[i,c(14:21)[good]]))
  pens[i] <- mean(as.numeric(z[i,c(6:13)[good]]))
}


inputs <- cbind(inputs,size,lengths,pens)
inputs <- as.matrix(inputs)
lm2 <- lm(allmeans[train] ~ inputs)
cor(lm2$fitted.values,allmeans[train])
plot(allmeans[train],lm2$fitted.values,xlab='Average Fraction Of True Rule Discovered',ylab='Fitted Fraction')

lm3 <- lm(completes ~ inputs)
plot(lm3$fitted.values,completes)
summary(lm3)

lm3 <- lm(fracCompletes ~ inputs)
plot(lm3$fitted.values,completes)
summary(lm3)

D <- as.data.frame(inputs)
lm1 = lm(allmeans[train] ~ inputs)
lm2 = lm(completes ~ inputs)
lm3 = lm(fracCompletes ~ inputs)


plot(allmeans[train],lm1$fitted.values,xlab='Average Fraction Of True Rule Discovered',ylab='Fitted Fraction')
plot(completes,lm2$fitted.values,xlab='Number of Complete Rules Discovered ',ylab='Fitted Cout')
plot(fracCompletes,lm3$fitted.values,xlab='Average Fraction Of True Rule Discovered',ylab='Fitted Fraction')


#### Examine the relationship between discovering true rules and model metrics & simulation inputs
bins = sort(unique(fracCompletes[-seq(10000)]))
big = seq(10000)
small = seq(20000)[-big]
meanROC = vector()
meanPR = vector()
features = vector()
for(i in 1:length(bins)){
  #hits = which(z$ROCAUC[seq(nrow(y))] <= bins[i] & z$ROCAUC[seq(nrow(y))] > bins[i-1])
  hits = which(fracCompletes == bins[i])
  #hits = intersect(hits,small)
  #meanComp[i] = mean(fracCompletes[hits])
  meanROC[i] = median(z$ROCAUC[hits])
  meanPR[i] = median(z$PRAUC[hits])
  features[i] = median(z$feats[hits])
  print(c(i,length(hits)))
  #points(bins[rep(i,length(hits))],z$ROCAUC[hits])
}

# FRACTION OF TRUE RULES RECOVERED
pdf(file='Fraction_of_rules_recovered.pdf')
par(mfrow=c(2,2))
plot(fracCompletes,z$ROCAUC,xlab='Fraction of True Rules Completely Recovered',ylab='Model ROC-AUC')
lines(bins,meanROC,col=2,lwd=2)
##legend("bottomright",#legend=c("Median ROC-AUC"),col=2,pch=16)


plot(fracCompletes,z$PRAUC,xlab='Fraction of True Rules Completely Recovered',ylab='Model PR-AUC')
lines(bins,meanPR,col=2,lwd=2)
##legend("bottomright",#legend=c("Median PR-AUC"),col=2,pch=16)

plot(fracCompletes,z$feats,ylab='Number of Features',xlab='Fraction of True Rules Completely Recovered',ylim=c(-5,1000))
lines(bins,features,col=2,lwd=2)
##legend("bottomright",#legend=c("Median Number of Features"),col=2,pch=16)

boxplot(fracCompletes[small],fracCompletes[big],ylab='Fraction of True Rules Completely Recovered',names=c('200','2000'),xlab='Sample Size')
dev.off()

mod=lm(fracCompletes ~ z$ROCAUC + z$PRAUC + z$feats + sampes)


library(lmvar)
X = cbind(fracCompletes,z$ROCAUC,z$PRAUC,z$feats,sampes)
fit = lm(fracCompletes ~ ., as.data.frame(X[,-1]), x = TRUE, y = TRUE)
cv.lm(fit)
#################################
# PROBABILITY OF RECOVERING AT LEAST 1 TRUE RULE
atleast1 = fracCompletes > 0
evalmod(scores=z$ROCAUC,labels = atleast1)
evalmod(scores=z$PRAUC,labels = atleast1)

bins = seq(20)/20
probs = vector()
probs[1]=0
fracRec = vector()
fracRec[1]=0


for(i in 2:length(bins)){
  hits = which(z$ROCAUC > bins[i-1] & z$ROCAUC <= bins[i])
  probs[i] = sum(atleast1[hits])/length(hits)
  fracRec[i] = mean(allmeans[hits])
}
plot(bins,probs)
plot(bins,fracRec)

disc = completes > 0
calibrate.plot(y = atleast1, p = z$ROCAUC)
mod = glm(atleast1 ~ z$ROCAUC + z$PRAUC + z$feats + sampes,family='binomial')
mod2 = randomForest(atleast1 ~.,  z$ROCAUC + z$PRAUC + z$feats + sampes)
library(popbio)
pdf(file='Probability_to_recover_minimum_1_rule.pdf')
par(mfrow=c(2,2))
logi.hist.plot(z$ROCAUC,atleast1,boxp=FALSE,type="hist",col="gray",xlab='Model ROC-AUC',main='Probability to Recover >= 1 True Rule')
logi.hist.plot(z$PRAUC,atleast1,boxp=FALSE,type="hist",col="gray",xlab='Model PR-AUC',main='Probability to Recover >= 1 True Rule')
logi.hist.plot(z$feats,atleast1,boxp=FALSE,type="hist",col="gray",xlab='# Features',main='Probability to Recover >= 1 True Rule')
sampes = c(rep(2000,10000),rep(200,10000))
logi.hist.plot(sampes,atleast1,boxp=FALSE,type="hist",col="gray",xlab='Sample Size',main='Probability to Recover >= 1 True Rule')
logi.hist.plot(mod$fitted.values,atleast1,boxp=FALSE,
               type="hist",col="gray",xlab='Fitted logistic regression probabilty ',main='Probability to Recover >= 1 True Rule')
dev.off()

### MEDIAN RANK OF RULES - WHEN RECOVERY DOES OCCUR
pdf('Median_Rank_when_discovered.pdf')
par(mfrow=c(2,2))
bad = which(bestrank > 5000) # abve i assigned an impossible rank to mark which simulations yielded nothing
plot(z$ROCAUC[-bad],log(bestrank[-bad]),ylab='Log Rank of First Rule',xlab='Model ROC-AUC')
plot(z$PRAUC[-bad],log(bestrank[-bad]),ylab='Log Rank of First Rule',xlab='Model PR-AUC')
plot(z$feats[-bad],log(bestrank[-bad]),ylab='Log Rank of First Rule',xlab='Number of Features')
boxplot(sampes[-bad],(bestrank[-bad]),ylab='Rank of First Rule',xlab='Number of Samples',names=c(200,2000))
dev.off()

rankMod = lm(log(bestrank[-bad]) ~ z$ROCAUC[-bad] + z$PRAUC[-bad] + z$feats[-bad] + sampes[-bad])

# AVERAGE PARTIAL RECOVERY OF RULES
pdf('Mean_Partial_Recovery.pdf')
par(mfrow=c(2,2))
plot(z$ROCAUC,allmeans[seq(20000)],xlab='Model ROC-AUC',ylab='Mean Partial Recovery')
plot(z$PRAUC,allmeans[seq(20000)],xlab='Model PR-AUC',ylab='Mean Partial Recovery')
plot(z$feats,allmeans[seq(20000)],xlab='Number of Features',ylab='Mean Partial Recovery')
boxplot(allmeans[10001:20000],allmeans[1:10000],xlab='Model ROC-AUC',ylab='Mean Partial Recovery',names=c(200,2000))
dev.off()

recMod = lm(allmeans[seq(20000)] ~z$ROCAUC + z$PRAUC + z$feats + sampes)


############# make a table of the correlations
correlationTable = matrix(0,ncol=4,nrow=4)
colnames(correlationTable) = c('ROC-AUC','PR-AUC','N Features','Sample Size')

correlationTable[1,1] = round(cor(z$ROCAUC,fracCompletes),digits=3)
correlationTable[1,2] = round(cor(z$PRAUC,fracCompletes),digits=3)
correlationTable[1,3] = round(cor(z$feats,fracCompletes),digits=3)
correlationTable[1,4] = round(cor(sampes,fracCompletes),digits = 3)

correlationTable[2,1] = round(cor(z$ROCAUC,means),digits=3)
correlationTable[2,2] = round(cor(z$PRAUC,means),digits=3)
correlationTable[2,3] = round(cor(z$feats,means),digits=3)
correlationTable[2,4] = round(cor(sampes,means),digits=3)

correlationTable[3,1] =  round(cor(z$ROCAUC[-bad],-bestrank[-bad]),digits=3)
correlationTable[3,2] =  round(cor(z$PRAUC[-bad],-bestrank[-bad]),digits=3)
correlationTable[3,3] =  round(cor(z$feats[-bad],-bestrank[-bad]),digits=3)
correlationTable[3,4] =  round(cor(sampes[-bad],-bestrank[-bad]),digits=3)

write.csv(correlationTable,file='CorrelationTable.csv')

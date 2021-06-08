

# Data for the mushroom dataset is found here: https://archive.ics.uci.edu/ml/machine-learning-databases/mushroom/

ruleMat = matrix(0,ncol=ncol(X),nrow=nrow(X))
j = seq(nrow(X))
matchMat = matrix(0,ncol=2,nrow=nrow(ruleMat))
best = 0
for(i in 1:length(j)){
  rules = JigSaw(RF = mod,data = X,traindata = X,labels = predicted,trainlabels = predicted,target = j[i],regularization = .005,categorical = categorical)
  ruleMat[j[i],rules[[1]]] = rules[[2]][1,]
  
  matches = seq(nrow(X))
  rule = rules[[1]]
  bounds  = rules[[2]]
  for(k in 1:length(rule)){
    matches = intersect(matches,which(X[,rule[k]] <= bounds[1,k] & X[,rule[k]] >= bounds[2,k]))
  }
  tab = table(labs[matches])
  matchMat[i,]= tab
  if(i %% 100 == 0){
    print(i)
  }
  if(tab[2] > best & tab[1] == 0){
    best = tab[2]
    print(best)
  }
}

vals = vector()
for(i in 1:nrow(ruleMat)){
  vals = union(vals,which(ruleMat[i,] != '0'))
}


write.csv(ruleMat,file='fungus_jigsaw_rules.csv')
#######################################
poisonlist = which(matchMat[,1]==0)
ord = rev(order(matchMat[poisonlist,1]))
poisonlist = poisonlist[ord]
count = 0
best = 0
lengths = vector()
vars = vector()
used = vector()
utilized = list()
vals = utilized
while(length(poisonlist) > 0){
  count = count + 1
  used = c(used,poisonlist[1])
  rules = JigSaw(RF = mod,data = X,traindata = X,labels = predicted,trainlabels = predicted,target = poisonlist[1],regularization = .005,categorical = categorical)
  ruleMat[poisonlist[1],rules[[1]]] = rules[[2]][1,]
  lengths = c(lengths,length(rules[[1]]))
  vars = union(vars,rules[[1]])
  
  matches = seq(nrow(X))
  rule = rules[[1]]
  utilized[[count]] = c(rule,matchMat[poisonlist[1],])
  
  bounds  = rules[[2]]
  vals[[count]]=rbind(rule,bounds[1,])
  for(k in 1:length(rule)){
    matches = intersect(matches,which(X[,rule[k]] <= bounds[1,k] & X[,rule[k]] >= bounds[2,k]))
  }
  tab = table(labs[matches])
  matchMat[poisonlist[1],]= tab
  if(count %% 100 == 0){
    print(count)
  }
  if(tab[1] > best & tab[2] == 0){
    best = tab[1]
    print(best)
  }
  poisonlist = poisonlist[-which(poisonlist %in% matches)]
}
count

#############
uRules = (unique(cbind(ruleMat,labs)))

rulelens = vector()
for(i in 1:nrow(uRules)){
  rulelens[i] = length(which(uRules[i,1:22]!='0'))
}

############ convert to 1s and 0s
tots = 0
for(i in 1:ncol(X)){
  tots = tots + length(table(X[,i]))
}

########## convert the features to binary predictors
binMat = matrix(0,ncol=tots,nrow=nrow(X))
binMat = as.data.frame(binMat)
miss = which(X[,11]=='?')
X[miss,11] = 'm'
featNmes=names(X)
matnames = vector()
origIndex = vector()
origCat = vector()
count = 0
for(i in 1:ncol(X)){
  categories = unique(X[,i])
  for(j in 1:length(categories)){
    count = count + 1
    hits = which(X[,i] == categories[j])
    binMat[hits,count] = 1
    matnames[count] = paste0(featNmes[i],categories[j],"_",i)
    origIndex[count] = i
    origCat[count] = categories[j]
    if(count == 56){
      print(i)
    }
  }
}

# convert the features to discretized features
discMat = matrix(0,ncol=ncol(X),nrow=nrow(X))
origCats = list()
for(i in 1:ncol(X)){
  categories = unique(X[,i])
  origCats[[i]] = categories
  for(j in 1:length(categories)){
    hits = which(X[,i]==categories[j])
    discMat[hits,i]=j
  }
}

labs2 = rep(0,length(labs))
# To 
labs2[which(labs=='p')]=1
labs2[which(labs=='e')]=0
labs2 = as.factor(labs2)
set.seed(10)
mod2 = randomForest(labs2 ~., binMat,keep.inbag=T,ntree=50)


source('Functions.r');h = BowSaw_ExistingSample(genes = binMat,mod = mod2,minmatch = 5,choice = 0)
tab = h[[4]][[2]]

prec = tab[,2]/tab[,1]
plot(tab[,2],prec)



j = which(prec == 1)


goodRules = unique(h[[3]])[j,]
goodRules = goodRules[-1,]


precG = vector()
caughtG = vector()
for(i in 1:nrow(goodRules)){
  matches = seq(nrow(binMat))
  rule = which(goodRules[i,]>0)
  vals = goodRules[i,rule]
  for(k in 1:length(rule)){
    matches = intersect(matches,which(binMat[,rule[k]] == vals[k]))
  }
  precG[i] = length(which(labs2[matches]==1))/length(matches)
  caughtG[i] = length(matches)
}

totrules = 0
for(i in 1:nrow(goodRules)){
  rule = which(goodRules[i,]>0)
  totrules = totrules + 2^length(rule) - 1
}

rulinos = matrix(0,ncol=ncol(goodRules),nrow=totrules)
valueinos = rulinos

count = 0
for(i in 1:nrow(goodRules)){
  
  rule = which(goodRules[i,]>0)
  values = goodRules[i,rule]
  combos = 2^length(rule)
  
  size = length(rule)
  if(size > 1){
    for(k in 1:size){
      news = combn(rule,k)
      for(m in 1:ncol(news)){
        count = count + 1
        rulinos[count,news[,m]]=1
        temp = which(rule %in% news[,1])
        valueinos[count,rule[temp]]=values[temp]
      }
    }
  }
  if(size==1){
    count = count + 1
    rulinos[count,rule]=1
  }
  
  
}

rulinos = rulinos[1:count,]
rulinos = unique(rulinos)

edibles = which(labs2==0)
found = list()
precisions = vector()
caught = vector()
foundMat = matrix(0,ncol=nrow(binMat),nrow=nrow(rulinos))
#foundMat[,which(labs2==1)]=2
for(i in 1:nrow(rulinos)){
  matches = seq(nrow(binMat))
  rule = which(valueinos[i,]>0)
  values = valueinos[i,rule]
  
  for(k in 1:length(rule)){
    matches = intersect(matches,which(binMat[,rule[k]] == values[k]))
  }
  found[[i]] = matches
  precisions[i] = length(which(labs2[matches]==0))/length(matches)
  caught[i] = length(matches)
  foundMat[i,matches]=1
  #edibles = edibles[-which(edibles %in% matches)]
  print(length(matches))
}

plot(caught,precisions)

all = vector()
good = which(precisions==1)
for(i in 1:length(good)){
  all = union(all,found[[good[i]]])
}
table(labs2[all])
table(labs2[-all])


###########
finders = foundMat[good,]
seekers = finders
used = vector()
added = vector()
applied = vector()
discovered = vector()


tops = 1
while(max(tops)>0){
  tops = rowSums(seekers)
  current = which.max(tops)
  if(max(tops)>0){
    hits = which(seekers[current,]==1)
    hits2 = which(finders[current,]==1)
    seekers[,hits]=0
    used = c(used,current)
    added = c(added,length(hits))
    applied = c(applied,length(hits2))
    discovered = union(discovered,hits)
  }

}
table(labs2[discovered])
table(labs2[-discovered])


## IF using binarized matrix
for(i in 1:length(used)){
  rule = which(valueinos[good[used[i]],]!=0)
  print(c(origIndex[rule],origCat[rule]))
}
applied
added

# if using categorical matrix
for(i in 1:length(used)){
  rules = which(valueinos[good[used[i]],]>0)
  vals = valueinos[good[used[i]],rules]
  
  cats = vector()
  for(k in 1:length(rules)){
    cats = c(cats,origCats[[rules[k]]][vals[k]])
  }
  print(c(rules,cats))
}
applied
added

plot(cumsum(added),ylab='Cumulative Finds')



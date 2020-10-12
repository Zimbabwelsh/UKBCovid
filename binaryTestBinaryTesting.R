library(tidyverse)

exposureType="binary"
partNum=192
numParts=200
traitofinterest = read.csv('C:/00RESEARCH/repo/UKB/testData/cneg33352.csv')
traitofinterestname = "cneg"
userId = "eid"

phenos = read.table(paste0('C:/00RESEARCH/repo/UKB/testData/data-', exposureType, '-', partNum, '-', numParts, '.txt'), sep=',', header=1, comment.char="")
confounders = data.frame(phenos$userID)


if (ncol(phenos)==1) {
  return(NULL)
}


## allow for estimation with no confounders
if (ncol(confounders)>=2){
  confNames = colnames(confounders)
  confNames = confNames[-which(confNames==userId)]
  
  ## merge datasets
  data = merge(traitofinterest, confounders, by=userId)
  data = merge(data, phenos, by.x=userId, by.y='userID')
  confs = data[,confNames]
} else{
  data=merge(traitofinterest, phenos, by.x=userId, by.y='userID')
  confs=confounders
}


## generate master list of phenos to loop over

phenoNames = colnames(phenos)
phenoNames = phenoNames[-which(phenoNames=='userID')]

for (i in 1:length(phenoNames)) {
  
  varName = phenoNames[i]
  
  print(varName)
  
  #sink(resLogFile, append=TRUE)
  #cat(paste(varName, ' || ', exposureType,' || ', sep=''))
  #sink()
  
  varType=''
  
  pheno = data[,varName]
  phenoFactor = factor(data[,traitofinterestname])
  
  ## if else catch loop for categorical data to specify reference categories and QC for n(exposure/not)
  
  if (exposureType != "cont"){
    pheno =as.factor(pheno)
    phenoF = chooseReferenceCategory(pheno)
    facLevels = levels(phenoF)
    
    cats = table(phenoF)
    idxMax = cats[which.max(cats)]
    idxMin = cats[which.min(cats)]
    numNotNA = length(which(!is.na(phenoF)))
    
    # if (idxMax<100) {
    #   print('<100')
    #   sink(resLogFile, append=TRUE)
    #   cat(paste('<100 (',idxMax,'/',idxMin,' || SKIP ',sep=''))
    #   sink()
    # }
    # else if (numNotNA<500) {
    #   print('<500 in total')
    #   sink(resLogFile, append=TRUE)
    #   cat('<500 in total (', numNotNA, ') || SKIP ',sep='')
    #   sink()
    # }
    # 
    # else {
    #   sink(resLogFile, append=TRUE)
    #   cat('testing || ')
    #   sink()
    # }
    
    ## specify exp is factor variable with most common category as reference
    
    exp = phenoF
    
    
    myBinaryRegression(exp, confs, phenoFactor, partNum, numParts, resDir, varName, exposureType)
    
    # save result to file
    # sink(resLogFile, append=TRUE)
    # cat(paste("SUCCESS results-", exposureType))
    # sink()
    # 
  }
  else {
    
    ## if non-categorical, exposure is just "pheno"
    
    exp = pheno
    
    myBinaryRegression(exp, confs, phenoFactor, partNum, numParts, resDir, varName, exposureType)
    
    # save result to file
    # sink(resLogFile, append=TRUE)
    # cat(paste("SUCCESS results-", exposureType))
    # sink()
  }
  
  # sink(resLogFile, append=TRUE)
  # cat("\n")
  # sink()
  # 

}

## regression function for flipped exp/phenoFactor

myBinaryRegression <- function(exp, confs, phenoFactor, partNum, numParts, resDir, varName, exposureType) {
  
  print(unique(phenoFactor))
  
  varType=exposureType
  
  ## binary logistic regression
  ## generalise to no confounders case
  if (ncol(confs)<2){
    mylogit <- glm(phenoFactor ~ exp, family="binomial")
  }
  
  else {
    mylogit <- glm(phenoFactor ~ exp + ., data=confs, family="binomial")
  }  
  
  sumx = coef(summary(mylogit))
  
  ## create if else loop for continuous, binary and un/ordered categotical
  ## for continuous generate N==numNotNA, exposure effects are single estimate
  
  if (exposureType == "cont"){
    pvalue = sumx['exp','Pr(>|z|)']
    beta = sumx["exp","Estimate"]
    cis = confint(mylogit, "exp", level=0.95)
    lower = cis["2.5 %"]
    upper = cis["97.5 %"]
    numNotNA = length(which(!is.na(exp)))
    
    
    # sink(resLogFile, append=TRUE)
    # cat(paste(varName, varType,  beta, lower, upper, pvalue, sep=","))
    # cat(' || ')
    # cat(paste(resDir,"results-",exposureType, "-", partNum, "-", numParts,".txt", sep=""))
    # sink()
    write(paste(varName, varType, numNotNA, beta, lower, upper, pvalue, sep=","), file=paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""), append=TRUE)
  }                                                                              
  
  ## for binary generate N = n(true)/n(false) (numNotNA), exposure still single estimates
  
  else if (exposureType == "binary") {
    pvalue = sumx['exp1','Pr(>|z|)']
    beta = sumx["exp1","Estimate"]
    cis = confint(mylogit, "exp1", level=0.95)
    lower = cis["2.5 %"]
    upper = cis["97.5 %"]
    facLevels = levels(exp)
    idxTrue = length(which(exp==facLevels[1]))
    idxFalse = length(which(exp==facLevels[2]))
    numNotNA = length(which(!is.na(exp)))
    
    # sink(resLogFile, append=TRUE)
    # cat(paste(varName, varType, paste(idxTrue,"/",idxFalse,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","))
    # cat(' || ')
    # cat(paste(resDir,"results-", exposureType,"-", partNum, "-", numParts,".txt", sep=""))
    # sink()
    write(paste(varName, varType, paste(idxTrue,"/",idxFalse,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","), file=paste("C:/00RESEARCH/repo/UKB/UKBCovid/results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""), append=TRUE)
  }
  
  ## for un/ordered categorical exposure - store multiple exposure effects, N = most common and numNotNA, loop over (facLevels-1)
  
  else {
    reference = levels(exp)[1]
    facLevels = levels(exp)
    cis       = confint(mylogit, level=0.95)
    cats      = table(exp)
    idxMax    = cats[which.max(cats)]
    numNotNA  = length(which(!is.na(exp)))
    
    for (i in (1:max(facLevels))){
      pvalue  = sumx[i+1,4]
      beta    = sumx[i+1,1]
      lower   = cis[i+1,1]
      upper   = cis[i+1,2]
      
      # sink(resLogFile, append=TRUE)
      # cat(paste(varName,i, "-", reference, varType, paste(idxTrue,"/",idxFalse,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","))
      # cat(' || ')
      # cat(paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""))
      # sink()
      write(paste(paste(varName, i,"-",reference, sep=""), varType, paste(maxFreq,"/",numNotNA,sep=""), beta, lower, upper, pvalue, sep=","), file=paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""), append=TRUE)
    }
  }
}

## chooseReferenceCategory taken from PHESANT "testCategoricalUnordered.r" 

chooseReferenceCategory <- function(pheno) {
  
  uniqVar = unique(na.omit(pheno));
  catFact = factor(pheno)
  
  maxFreq=0;
  maxFreqVar = "";
  for (u in uniqVar) {
    withValIdx = which(pheno==u)
    numWithVal = length(withValIdx);
    if (numWithVal>maxFreq) {
      maxFreq = numWithVal;
      maxFreqVar = u;
    }
  }
  
  cat("reference: ", maxFreqVar,"=",maxFreq, " || ", sep="");
  
  ## choose reference (category with largest frequency)
  phenoFactor <- relevel(pheno, ref = paste("",maxFreqVar,sep=""))
  
  return(phenoFactor);
}

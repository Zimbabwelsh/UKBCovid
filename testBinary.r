testBinary <- function(resDir, partNum, numParts, confounders, traitofinterest, traitofinterestname, userId, phenoDir, exposureType) {

  ## create LogFile
  # resLogFile = paste(resDir, "results-log-", exposureType, "-", partNum, "-", numParts, ".txt",sep="")
  # sink(resLogFile)
  # sink()
  
  ## generate empty results file
  # write("varName,varType,n,beta,lower,upper,pvalue", file=paste(opt$resDir,"results-", exposureType, "-", opt$partIdx, "-", opt$numParts, ".txt",sep=""), append=FALSE)
  
  
  resLogFile = paste(resDir, "results-log-", exposureType, "-", partNum, "-", numParts, ".txt",sep="")
  
  ## load phenotype data
  phenos = read.table(paste(phenoDir, '/data-', exposureType,"-", partNum, '-', numParts, '.txt',sep=''), sep=',', header=1, comment.char="")
	
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
	} 
  else {
	  data=merge(traitofinterest, phenos, by.x=userId, by.y='userID')
	  confs=confounders
	}
	
	
	## generate master list of phenos to loop over
	
	phenoNames = colnames(phenos)
	phenoNames = phenoNames[-which(phenoNames=='userID')]
	
	for (i in 1:length(phenoNames)) {
	  
	  varName = phenoNames[i]
	  
	  print(varName)
	  
	  sink(resLogFile, append=TRUE)
	  cat(paste(varName, ' || ', exposureType,' || ', sep=''))
	  sink()
	  
	  varType=''
	  
	  ## generate variables for inclusion in regression
	  
	  # pheno is exposure phenotype
	  exp = data[,varName]
	  
	  # caseness is binary outcome for subsample
	  caseness = factor(data[,traitofinterestname])
	  
	  ## if else catch loop for categorical data to specify reference categories and QC for n(exposure/not)
	  ## Keep consistent - continuous, then binary, then *else*
	  
	  # If non-categorical, exposure is just "pheno"
	  
    	  if (exposureType == "cont"){
    	    
    	    numNotNA = length(which(!is.na(exp)))
    	    
    	    if (numNotNA<500){
    	      sink(resLogFile, append=TRUE)
    	      cat('<500 in total (', numNotNA, ') || SKIP', sep='')
    	      sink()
    	    }
    	    else{
    	      myBinaryRegression(caseness, confs, exp, partNum, numParts, resDir, varName, exposureType)
    	    
    	    # save result to file
    	    sink(resLogFile, append=TRUE)
    	    cat(paste("SUCCESS results-", exposureType))
    	    sink()
    	    }
    	  }
	  
	  # For binary exposure, read pheno as factor and exclude if cases/controls less than 10 or total less than 500.
	  
	      else if (exposureType == "binary"){
	    
  	    exp =as.factor(exp)
  	    facLevels = levels(exp)
  	    idxTrue = length(which(exp==facLevels[1]))
  	    idxFalse = length(which(exp==facLevels[2]))
  	    numNotNA = length(which(!is.na(exp)))
      	    
      	    if (idxTrue<10 || idxFalse<10) {
      	      print('<10')
      	      sink(resLogFile, append=TRUE)
      	      cat(paste('<10 (',idxTrue,'/',idxFalse,' || SKIP ',sep=''))
      	      sink()
      	    }
      	    else if (numNotNA<500) {	
      	      print('<500 in total')
      	      sink(resLogFile, append=TRUE)
      	      cat('<500 in total (', numNotNA, ') || SKIP ',sep='')
      	      sink()
      	    }
      	    else {
      	      
      	      sink(resLogFile, append=TRUE)
      	      cat('testing || ')
      	      sink()
      	      
      	      myBinaryRegression(caseness, confs, exp, partNum, numParts, resDir, varName, exposureType)
      	      
      	      ## save result to file
      	      sink(resLogFile, append=TRUE)
      	      cat("SUCCESS results-binary")
      	      sink()
      	    }
	      }
	      else {
	          exp =as.factor(exp)
	          exp = chooseReferenceCategory(exp)
	          facLevels = levels(exp)
	          
	          cats = table(exp)
	          idxMax = cats[which.max(cats)]
	          idxMin = cats[which.min(cats)]
	          numNotNA = length(which(!is.na(exp)))
	          
	          if (idxMax<100) {
	            print('<100')
	            sink(resLogFile, append=TRUE)
	            cat(paste('<100 (',idxMax,'/',idxMin,' || SKIP ',sep=''))
	            sink()
	          }
	          else if (numNotNA<500) {
	            print('<500 in total')
	            sink(resLogFile, append=TRUE)
	            cat('<500 in total (', numNotNA, ') || SKIP ',sep='')
	            sink()
	          }
	          
	          else {
	            sink(resLogFile, append=TRUE)
	            cat('testing || ')
	            
	            #include variable generation in the else loop - so only generate analysis where n is reasonable.
	            
	            myBinaryRegression(caseness, confs, exp, partNum, numParts, resDir, varName, exposureType)
	            
	            # save result to file
	            sink(resLogFile, append=TRUE)
	            cat(paste("SUCCESS results-", exposureType))
	            sink()
	            
	            sink()
	          }
	      }
	          
  	  
    	  sink(resLogFile, append=TRUE)
    	  cat("\n")
    	  sink()
  	  }
}
	
	## regression function for flipped exp/phenoFactor
	
	myBinaryRegression <- function(caseness, confs, exp, partNum, numParts, resDir, varName, exposureType) {
	  
	  print(unique(caseness))
	  
	  resLogFile = paste(resDir, "results-log-", exposureType, "-", partNum, "-", numParts, ".txt",sep="")
	  
	  varType=exposureType
	  
	  ## binary logistic regression
	    ## generalise to no confounders case
	  if (ncol(confs)<2){
	    mylogit <- glm(caseness ~ exp, family="binomial")
	  }
	  
	  else {
	    mylogit <- glm(caseness ~ exp + ., data=confs, family="binomial")
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
	    
	    
	    sink(resLogFile, append=TRUE)
	    cat(paste(varName, varType,  beta, lower, upper, pvalue, sep=","))
	    cat(' || ')
	    cat(paste(resDir,"results-",exposureType, "-", partNum, "-", numParts,".txt", sep=""))
	    sink()
	    
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
	    
	    sink(resLogFile, append=TRUE)
	    cat(paste(varName, varType, paste(idxTrue,"/",idxFalse,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","))
	    cat(' || ')
	    cat(paste(resDir,"results-", exposureType,"-", partNum, "-", numParts,".txt", sep=""))
	    sink()
	    
	    write(paste(varName, varType, paste(idxTrue,"/",idxFalse,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","), file=paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""), append=TRUE)
	  }
	  
	  ## for un/ordered categorical exposure - store multiple exposure effects, N = most common and numNotNA, loop over 1:max(facLevels-1)
	  ## append row for each category.
	  
	  else {
	    reference = levels(exp)[1]
	    facLevels = levels(exp)
	    cis       = confint(mylogit, level=0.95)
	    cats      = table(exp)
	    idxMax    = cats[which.max(cats)]
	    idxMin    = cats[which.min(cats)]
	    numNotNA  = length(which(!is.na(exp)))
	    
	    for (i in 0:(as.numeric(max(facLevels))-1)){
	      pvalue  = sumx[i+2,4]
	      beta    = sumx[i+2,1]
	      lower   = cis[i+2,1]
	      upper   = cis[i+2,2]
	     
	      sink(resLogFile, append=TRUE)
	      cat(paste(varName, "-", i, "-", varType, paste(idxMax,"/",idxMin,"(",numNotNA,")",sep=""), beta, lower, upper, pvalue, sep=","))
	      cat(' || ')
	      cat(paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""))
	      sink()
	    
	      write(paste(paste(varName, "-", i, sep=""), varType, paste(idxMax,"/", idxMin, "(", numNotNA,")", sep=""), beta, lower, upper, pvalue, sep=","), file=paste(resDir,"results-", exposureType, "-", partNum, "-", numParts,".txt", sep=""), append=TRUE)
	    }
	  }
}
	
	## Function to select most common category as reference for unordered and ordered categorical variables:
	## chooseReferenceCategory taken from PHESANT/PHESANT-from-saved/ "testCategoricalUnordered.r" 
	
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
	


library("optparse")

option_list = list(
  make_option(c("-e", "--traitofinterest"), type="character", default=NULL, help="traitofinterest option should specify trait of interest variable name", metavar="character"),
  make_option(c("-g", "--traitofinterestfile"), type="character", default=NULL, help="Trait of interest dataset file name", metavar="character"),
  make_option(c("-r", "--resDir"), type="character", default=NULL, help="resDir option should specify directory where results files should be stored", metavar="character"),
  make_option(c("-u", "--userId"), type="character", default="userId", help="userId option should specify user ID column in trait of interest and phenotype files [default= %default]", metavar="character"),
  make_option(c("-a", "--partIdx"), type="integer", default=NULL, help="Part index of phenotype (used to parellise)"),
  make_option(c("-b", "--numParts"), type="integer", default=NULL, help="Number of phenotype parts (used to parellise)"),
  make_option(c("-c", "--confounderfile"), type="character", default=NULL, help="Confounder file name", metavar="character"),
  make_option(c("-p", "--phenoDir"), type="character", default=NULL, help="Phenotype directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

source('testBinary.r')

## load confounders
confounders = read.table(opt$confounderfile, sep=',', header=1)

## load trait of interest
exposure = read.table(opt$traitofinterestfile, sep=',', header=1)

# move logfile and results creation to testBinary.r

## create logFile
# resLogFile = paste(opt$resDir, "results-log-", opt$partIdx, "-", opt$numParts, ".txt",sep="")
# sink(resLogFile)
# sink()


## create empty results file

# write("varName,varType,n,beta,lower,upper,pvalue", file=paste(opt$resDir,"results-logistic-binary-", opt$partIdx, "-", opt$numParts, ".txt",sep=""), append=FALSE)

## test each type of outcome

testBinary(opt$resDir, opt$partIdx, opt$numParts, confounders, exposure, opt$traitofinterest, opt$userId, opt$phenoDir, "cont")
testBinary(opt$resDir, opt$partIdx, opt$numParts, confounders, exposure, opt$traitofinterest, opt$userId, opt$phenoDir, "catunord")
testBinary(opt$resDir, opt$partIdx, opt$numParts, confounders, exposure, opt$traitofinterest, opt$userId, opt$phenoDir, "catord")
testBinary(opt$resDir, opt$partIdx, opt$numParts, confounders, exposure, opt$traitofinterest, opt$userId, opt$phenoDir, "binary")





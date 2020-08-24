#!/bin/bash

#request resources:
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -o output
#PBS -e errors

cd $PBS_O_WORKDIR

module load languages/R-3.6.2-gcc9.1.0

date

dataDir="${HOME}/data/UKB/15825/2019-05-02/"
codeDir="${HOME}/data/UKB/15825/2019-05-02/PHESANT/WAS/"
varListDir="${HOME}/data/UKB/15825/2019-05-02/PHESANT/variable-info/"

outcomeFile="${dataDir}data/data.33352.csv"
varListFile="${varListDir}outcome-info.tsv"
dcFile="${varListDir}data-coding-ordinal-info.txt"
resDir="${dataDir}derived"

# start and end index of phenotypes
pIdx=1
np=10

cd $codeDir
Rscript ${codeDir}phenomeScan.r --partIdx=$pIdx --numParts=$np --phenofile=${outcomeFile} --variablelistfile=${varListFile} --datacodingfile=${dcFile} --resDir=${resDir} --userId="eid" --save



date

#!/bin/bash
#PBS -l walltime=360:00:00,nodes=1:ppn=1
#PBS -o output
#PBS -e errors
#PBS -t 1-100
#---------------------------------------------

datadir="${HOME}/data/UKB/15825/2019-05-02/derived"
origdir="${HOME}/data/UKB/15825/2019-05-02/data/"
head -n 1 ${origdir}data.33352.csv | sed 's/,"/,"x/g' | sed 's/-/_/g' | sed 's/\./_/g' > ${datadir}data.33352-phesant_header.csv
awk '(NR>1) {print $0}' ${origdir}data.33352.csv >> ${datadir}data.33352-phesant_header.csv

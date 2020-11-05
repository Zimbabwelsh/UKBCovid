library(tidyverse)
library(lubridate)
library(data.table)

## Remove eids for deaths prior to January 2020
deaths <- read.table("data/death-20200916.txt", h=TRUE)
deaths$date_of_death <-  dmy(deaths$date_of_death)
deaths <- deaths %>% filter(date_of_death<"2020-01-01")


## Read in data for assessment centre - remove testing centres outside of England due to PHE data.
df <- fread("data/data.43017.phesant.tab", select=c('eid', 'x54_0_0'), he=T)
df <- df %>% mutate( england = case_when(x54_0_0 != 11003 &
                                           x54_0_0 != 11005 &
                                           x54_0_0 != 11004 &
                                           x54_0_0 != 11022 &
                                           x54_0_0 != 11023 ~ 1,
                                         (TRUE~0)))
df <- df[england==1]

## Read in withdrawals and omit from data
wd <- fread("data/withdrawals/meta.withdrawn.20200820.csv")
deathswds <- as.data.frame(append(wd$V1, deaths$eid))
colnames(deathswds) <- "eid"

## Create and write out master list of eligible EIDs (alive by Jan 2020, and English assessment centre, not withdrawn)
list <- subset(df, !(eid %in% deathswds$eid))
write.table(list$eid, file="data/phewasEids.txt",col.names = "eid", row.names=FALSE)

## Read in COVID test results
covid <- fread("data/covid19_results-20200909.txt", h=TRUE)

# Generate mean test result - treating "positive" as "ever testing positive" versus "never testing positive"
covid <- covid %>% group_by(eid) %>% mutate(postests=mean(result))
covid$cneg <- as.numeric(covid$postests==0)
covid$cpos <- as.numeric(covid$postests!=0)

postest <- setDT(unique(covid %>% select('eid','cpos')),key=NULL, check.names=FALSE)
negtest <- setDT(unique(covid %>% select('eid','cneg')), key=NULL, check.names=FALSE)

# Combine testing data with wider UKB cohort filtered by eligibility
df2<-merge(df, negtest, all=TRUE, by="eid")
df2$cneg[is.na(df2$cneg)]<-0
df2 <- subset(df2, (eid %in% list$eid))

df3 <- merge(df, postest, all=TRUE, by="eid")
df3$cpos[is.na(df3$cpos)]<-0
df3 <- subset(df3, (eid %in% list$eid))

## Write out traitofinterest files for PHESANT input
write.table(df2 %>% select(c("eid", "cneg")), file="data/cneg.csv", col.names=c("eid", "cneg"), row.names=FALSE)
write.table(df3 %>% select(c("eid", "cpos")), file="data/cpos.csv", col.names=c("eid", "cpos"), row.names=FALSE)

## Read in age and sex data for generating confounder files filtered by eligibility. 
confs <- fread("data/data.43017.phesant.tab", select=c('eid','x31_0_0', 'x34_0_0'), he=T)
confs <- subset(confs,(eid %in% list$eid))
write.table(confs, file="data/confs.csv", row.names=FALSE)



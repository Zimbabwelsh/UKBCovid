library(data.table)
library(dplyr)
library(parallel)

# read in covid testing
covid <- fread("data/covid_19/covid19_result_2020_05_26.txt", he=T)
# create phenotypes of positive and negative test
covid <- covid %>% group_by(eid) %>% mutate(postests=mean(result))
covid$testsneg <- as.numeric(covid$postests==0)
covid$testspos <- as.numeric(covid$postests!=0)
postest <- setDT(unique(covid %>% select('eid','testspos')),keep.rownames=TRUE, key=NULL, check.names=FALSE)
postest <- postest[testspos==1]
negtest <- setDT(unique(covid %>% select('eid','testsneg')),keep.rownames=TRUE, key=NULL, check.names=FALSE)
negtest <- negtest[testsneg==1]


# create subset of English patients who are alive at start of covid testing.
# testing data only available for England.
df <- fread("data/data.33352.csv", select=c('eid','54-0.0', '40007-0.0', he=T))
# generate vars for "english" and "deceased"
df <- df %>% mutate(deceased= case_when(Deathage0==TRUE~1, 
                                        is.na(Deathage0)~0,
                                        Deathage1==TRUE~1,
                                        Deathage2==TRUE~1),
                    england = case_when(ukbass != 11003 &
                                          ukbass != 11005 &
                                          ukbass != 11004 &
                                          ukbass != 11022 &
                                          ukbass != 11023 ~ 1,
                                        (TRUE~0))
                    )
# omit non-english and deceased
df <- df[england==1 & deceased==0]

# Generate case definitions
df$cpos <- as.numeric(df$eid %in% postest$eid)
df$cneg <- as.numeric(df$eid %in% negtest$eid)

df2 <- fread("data/data.33352.csv", he=T)
df3 <- left_join(df[,c("eid", "cpos")], df2, by ='eid')
df4 <- left_join(df[,c("eid", "cneg")], df2, by ='eid')

write.csv(df3,"33352cpos.csv", row.names = TRUE)
write.csv(df4,"33352cneg.csv", row.names = TRUE)
            
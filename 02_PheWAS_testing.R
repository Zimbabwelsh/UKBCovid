library(data.table)
library(dplyr)
library(parallel)

# read in covid testing
covid <- fread("data/covid_19/covid19_result_2020_07_31.txt", he=T)

#create phenotypes of positive and negative test
covid <- covid %>% group_by(eid) %>% mutate(postests=mean(result))
covid$testsneg <- as.numeric(covid$postests==0)
covid$testspos <- as.numeric(covid$postests!=0)
postest <- setDT(unique(covid %>% select('eid','testspos')),keep.rownames=TRUE, key=NULL, check.names=FALSE)
postest <- postest[testspos==1]
negtest <- setDT(unique(covid %>% select('eid','testsneg')),keep.rownames=TRUE, key=NULL, check.names=FALSE)
negtest <- negtest[testsneg==1]


#Indicate source location for phenotypes from PHESANT
#Need a way to string together predictors across PHESANT .txt files
#OR need to manipulate all PHESANT output so each file is a single predictor. 
#Then can use Gib's code directly for a$discovery, assuming each column is titled $discovery

# NEED to make sure that the non-testing regions are omitted, generate list of EIDs from non-England,

df <- df[england==1]
df <- df[deceased==0]

# Also remove deceased and withdrawn consent (think maybe done already?)
# Or run Louise's PHESANT package directly from command line? Need to pre-process to omit results from outside England 
# and withdrawn/deceased participants

phenpath <- "data/derived/phesant_mod"
phens <- list.files(phenpath)


phenout <- mclapply(phens, function(x)
  {
  message(x)
  fn <- file.path(phenpath, x, "phen.txt")
  file.exists(fn)
  a <- fread(fn)
  a <- merge(negtest, a, by.x="eid", by.y="FID")
  form <- paste0("negtest~discovery")
  return(try(summary(glm(form, a, family="binomial"))$coefficients))
}, mc.cores=6)

#


# Remove non-testing regions


# Generate case definitions
df$covidpos <- as.numeric(df$eid %in% postest$eid)
df$covidneg <- as.numeric(df$eid %in% negtest$eid)

#read in UKB data
df <- fread("data/data.33352.csv", select=c('eid', '31-0.0', '34-0.0', '54-0.0', '93-0.0', '94-0.0', '189-0.0',
                                            '884-0.0', '894-0.0', '904-0.0', '914-0.0', '1558-0.0', '1568-0.0', '1578-0.0',
                                            '1588-0.0', '1598-0.0', '2443-0.0', '4079-0.0', '4080-0.0', '6138-0.0',
                                            '6150-0.0', '21000-0.0', '21001-0.0', '20116-0.0', '40007-0.0', '40007-1.0', '40007-2.0'), he=T)

#rename vars
df=df %>% rename(male = "31-0.0", yob = '34-0.0', ukbass= '54-0.0', mansysbp= '93-0.0', mandiabp='94-0.0', towns='189-0.0',
                 mpaweek='884-0.0',mpatime= '894-0.0',vpaweek= '904-0.0',vpatime= '914-0.0', alccons = '1558-0.0', redwine='1568-0.0',whtwine= '1578-0.0',
                 beer='1588-0.0',spirits= '1598-0.0',diabetes= '2443-0.0',autdiabp= '4079-0.0',autsysbp ='4080-0.0', qualif='6138-0.0',
                 vascprob='6150-0.0',ethnic= '21000-0.0', bmi='21001-0.0', smoking='20116-0.0', Deathage0='40007-0.0', Deathage1='40007-1.0', Deathage2='40007-2.0')

#remove missingness
df <- df[smoking<3 & smoking>=0]
df <- df[ethnic>=0]
df <- df[!is.na(qualif)]

#gen derived vars
df <- df %>% mutate(highbp = case_when((autsysbp<140 & autdiabp<90) ~ 0,
                                       (autsysbp>=140 & !is.na(autdiabp)) ~ 1,
                                       (autdiabp>=90 & !is.na(autsysbp)) ~ 1),
                    
                    
                    alc = case_when(alccons>=4 & alccons<=6 ~ 0,
                                    alccons>=1 & alccons<=3 ~ (whtwine+redwine)*2+beer+spirits),
                    
                    
                    hialc = case_when((alc>=14 & male ==0) ~ 2,
                                      (alc>=21 & male ==1) ~ 2,
                                      (alc>0 & alc<14 & male==0 ~ 1),
                                      (alc>0 & alc<21 & male==1 ~ 1),
                                      alc==0~0),
                    
                    active= case_when((vpatime*vpaweek >= 75 | mpatime*mpaweek >= 150)~0,
                                      (mpaweek==0 ~ 2),
                                      (TRUE~1)),
                    
                    obesity= case_when(bmi<25~0,
                                       bmi<=25 & bmi <30 ~1,
                                       bmi<=30 ~2),
                    
                    nwhite= case_when(ethnic==1001 | ethnic==1002 | ethnic == 1003 | ethnic==1 ~0,
                                      (TRUE~1)),
                    
                    deceased= case_when(Deathage0==TRUE~1, 
                                        is.na(Deathage0)~0,
                                        Deathage1==TRUE~1, Deathage2==TRUE~1
                    ),
                    
                    england = case_when(ukbass != 11003 &
                                          ukbass != 11005 &
                                          ukbass != 11004 &
                                          ukbass != 11022 &
                                          ukbass != 11023 ~ 1,
                                        (TRUE~0)),
                    
                    lifestyle = hialc + active + obesity + smoking
                    
) 

# Remove non-testing regions
df <- df[england==1]
df <- df[deceased==0]

# Generate case definitions
df$covidpos <- as.numeric(df$eid %in% postest$eid)
df$covidneg <- as.numeric(df$eid %in% negtest$eid)




####


nwhitepos <- summary(glm(covidpos ~ nwhite, df, family="binomial"))$coefficients[2,1:2]

nwhiteneg <- summary(glm(covidneg ~ nwhite, df, family="binomial"))$coefficients[2,1:2]

smokepos <- summary(glm(covidpos ~ smoking, df, family="binomial"))$coefficients[2,1:2]

smokeneg <- summary(glm(covidneg ~ smoking, df, family="binomial"))$coefficients[2,1:2]

townspos <- summary(glm(covidpos ~ towns, df, family = "binomial"))$coefficients[2,1:2]

townsneg <- summary(glm(covidneg ~ towns, df, family = "binomial"))$coefficients[2,1:2]

obespos <- summary(glm(covidpos ~ obesity, df, family = "binomial"))$coefficients[2,1:2]

obesneg <- summary(glm(covidneg ~ obesity, df, family = "binomial"))$coefficients[2,1:2]

data <- rbind(nwhitepos, nwhiteneg, smokepos, smokeneg, townspos, townsneg, obespos, obesneg)
data <- as.data.frame(data)
data$OR <- exp(data[,1])
data$LoCI <- exp(data[,1])-(1.96*(data[,2]))
data$HiCI <- exp(data[,1])+(1.96*(data[,2]))

write.csv(data,"posneg.csv", row.names = TRUE)


###
# Attempt to write function returning model ests for pos and neg tests

#varlist <- list('nwhite', 'smoking', 'towns', 'obesity', 'male', 'bmi')
#results <- list()

#results <- lapply(varlist, function(x) {
#  summary(glm(covidpos ~ varlist, df, family="binomial"))$coefficients[2,1:2]
#})

#results <- do.call(rbind, results)
#results


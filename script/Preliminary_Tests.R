## ---------------------------
## Author: Helene Cecilia
##
## Date Created: 2023-06-12

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(chron)
library(glmmTMB)
library(DHARMa)
library(effects)
library(MuMIn)
library(multcomp) # for glht
library(broom) # for tidy

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Helper functions ----
hms_to_decimal_day <-function(time_hms){
  h <- hours(time_hms)
  m <- minutes(time_hms)
  s <- seconds(time_hms)
  sum_minutes <- h*60 + m + s/60
  dec <- sum_minutes / (24*60)
  return(dec)
}

# Load data ----
df1 <- read.csv("../data/DENV_Cyno_data_engorgement.csv")
df1$NHP <- "Cyno"
df1$virus <- "Dengue virus"

df2 <- read.csv("../data/DENV_Squirrel_data_engorgement.csv")
df2$NHP <- "Squirrel"
df2$virus <- "Dengue virus"

df3 <- read.csv("../data/ZIKV_Squirrel_data_engorgement.csv")
df3$NHP <- "Squirrel"
df3$virus <- "Zika virus"
# removing controls as they're already in df2
df3 <- df3[df3$group != "control",]

df4 <- read.csv("../data/ZIKV_Cyno_data_engorgement.csv")
df4$NHP <- "Cyno"
df4$virus <- "Zika virus"

df <- rbind(df1,df2,df3,df4)

df$monkey_status <- df$virus
df$monkey_status[which(df$group == "control")] <- "Control"

df$group_agg <- interaction(df$NHP,df$monkey_status)
df$group_agg <- gsub(" ","_",df$group_agg) # no space for posthoc test to work
df$group_agg <- as.factor(df$group_agg)

feed_duration <- str_split(df$feed_duration,":")
vec_duration <- sapply(feed_duration,"[[",2)
df$int_duration <- as.integer(vec_duration)

df$is_day28 <- df$day == 28
df$is_day28 <- as.factor(df$is_day28)
df$time_of_day <- df$dec - df$day

# Test if day 28 different from other days ----
## Duration ----
### Cyno ----
test <- df[df$NHP == "Cyno",]

summary(test$int_duration[test$is_day28 == TRUE])
summary(test$int_duration[test$is_day28 == FALSE])

m1 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = poisson) 
# model convergence problem with nbinom1 (account for overdispersion), so stick with poisson
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 723.037 

m2 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
              data = test,
              family = poisson) 
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # signif
MuMIn::AICc(m2) # 724.308 

anova(m1,m2) # not signif

# m1 is the best we can do
plot(allEffects(m1, partial.residuals = T))
summary(m1) # no effect of day 28
car::Anova(m1, type = "II") # no effect of day 28

### Squirrel ----
test <- df[df$NHP == "Squirrel",]

summary(test$int_duration[test$is_day28 == TRUE])
summary(test$int_duration[test$is_day28 == FALSE])

m1 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = poisson) 
# model convergence problem with nbinom1 (account for overdispersion), so stick with poisson
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 698.1084

m2 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
              data = test,
              family = poisson) 
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # signif
MuMIn::AICc(m2) # 700.1834

anova(m1,m2) # not signif

# m1 is the best we can do
plot(allEffects(m1, partial.residuals = T))
summary(m1) # effect of day 28
car::Anova(m1, type = "II") # effect of day 28

# longer durations on day 28
res <- confint(m1)
exp(res) 

## Time of day ----
# Some summary stats
summary(df$time_of_day)
(0.6813*24*60)/60 # 16.3512 hours
0.3512*60 # 21.1 minutes
summary(df$time_of_day[df$is_day28 == FALSE])
(0.3618056*24*60)/60 # 8.683334 hours
0.683334*60 #41 minutes
(0.6038194*24*60)/60 # 14.49167 hours
0.49167*60 # 29.5 minutes
df_no28 <- df[df$is_day28 == FALSE,]
hist(df$time_of_day)

### Cyno  ----
test <- df[df$NHP == "Cyno",]
# Time of day is between 0 and 1, but 
# i) as it's a cycle, 0 and 1 are actually close in time 
# and ii) our values don't go near the extreme so we'll use linear models 

m1 <- glmmTMB(time_of_day ~ is_day28,
              data = test) 
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
MuMIn::AICc(m1) # -739.6126 

m2 <- glmmTMB(time_of_day ~ is_day28 + (1|ID),
              data = test)
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
MuMIn::AICc(m2) # -863.4266 

anova(m1,m2) # in favour of random effects

plot(allEffects(m2, partial.residuals = T))
summary(m2) # no effect of day 28
car::Anova(m2, type = "II") # no effect of day 28

### Squirrel ----
test <- df[df$NHP == "Squirrel",]

m1 <- glmmTMB(time_of_day ~ is_day28,
              data = test) 
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # -457.4701

m2 <- glmmTMB(time_of_day ~ is_day28 + (1|ID),
              data = test)
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # ok
MuMIn::AICc(m2) # -455.371

anova(m1,m2) # NOT in favour of random effects

plot(allEffects(m1, partial.residuals = T))
summary(m1) # effect of day 28
# later in the day on day 28
res <- confint(m1)

# At this point we confirmed the need to exclude day 28 from all analyses

# Effect of species and/or infection group -----
# to assess the need to include interaction terms
# Without day 28 

## Weight (after day 0) ----

m1 <- glmmTMB(weight ~ group_agg,
              data = df[df$day != 28,])
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # 4384.206

# select only the biologically relevant comparisons
tuk = glht(m1, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # effect of species (both) and infection in cyno
summary(tuk, test = adjusted("fdr")) # same

plot(allEffects(m1, partial.residuals = T))
res <- tidy(confint(tuk))
# lower weight in DENV-infected cyno than control

## Temperature -----

### Temperature at the time of feeding (after day 0) ----

# Without day 28 

df <- df[df$day != 28,]

m1 <- glmmTMB(temp_estim_feed ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
MuMIn::AICc(m1) # 350.8832

m2 <- glmmTMB(temp_estim_feed ~ group_agg + (1|ID) + (1|day),
              data = df)
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
MuMIn::AICc(m2) # 286.7474

anova(m1, m2) # in favour of random effects

# select only the biologically relevant comparisons
tuk = glht(m2, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Zika_virus - Squirrel.Zika_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Cyno.Zika_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0")))
summary(tuk, test = adjusted("none")) # species diff both signif + infected/control in cyno
summary(tuk, test = adjusted("fdr")) # species diff both signif + infected/control in cyno (close to 0.05)

plot(allEffects(m2, partial.residuals = T))
# higher temp in DENV-infected cynos than control, but CI contains 0
res <- tidy(confint(tuk))

### All temperature data (every 15 mins) ----
df1 <- read.csv("../data/Temperatures_Sylvatic_DENV-2_Cynomolgus_Macaques.csv")
df2 <- read.csv("../data/Temperatures_Sylvatic_DENV-2_Squirrel_Monkeys.csv")
df3 <- read.csv("../data/Temperatures_Sylvatic_ZIKV_Squirrel_Monkeys.csv")
df4 <- read.csv("../data/Temperatures_Sylvatic_ZIKV_Cynomolgus_Macaques.csv")

treat1 <- read.csv("../data/Hanley2024_Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                   dec = ".", sep  ="\t")
treat2 <- read.csv("../data/Hanley2024_Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat3 <- read.csv("../data/Hanley2024_Table_S4_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat4 <- read.csv("../data/Hanley2024_Table_S3_Sylvatic_ZIKV_Cynomolgus_Macaques.csv",
                   dec = ".", sep  ="\t")

treat1 <- unique(treat1[,c("ID","Final.Treatment")])
treat2 <- unique(treat2[,c("ID","Final.Treatment")])
treat3 <- unique(treat3[,c("ID","Final.Treatment")])
treat4 <- unique(treat4[,c("ID","Final.Treatment")])

df1 <- merge(df1,treat1, by = "ID")
df2 <- merge(df2,treat2, by = "ID")
df3 <- merge(df3,treat3, by = "ID")
df4 <- merge(df4,treat4, by = "ID")

# merge controls for same species experiments
cont1 <- df1[df1$Final.Treatment == "Control",]
cont2 <- df2[df2$Final.Treatment == "Control",]
cont3 <- df3[df3$Final.Treatment == "Control",]

df2 <- rbind(df2,cont3)
df3 <- rbind(df3,cont2)
df4 <- rbind(df4,cont1)

df1$chron_time <- chron(times. = as.character(df1$time),
                        format = "h:m:s")
df2$chron_time <- chron(times. = as.character(df2$time),
                        format = "h:m:s")
df3$chron_time <- chron(times. = as.character(df3$time),
                        format = "h:m:s")
df4$chron_time <- chron(times. = as.character(df4$time),
                        format = "h:m:s")

df1$group <- df1$Final.Treatment
df2$group <- df2$Final.Treatment
df3$group <- df3$Final.Treatment
df4$group <- df4$Final.Treatment

df1$group[df1$group != "Control"] <- "Infected"
df2$group[df2$group != "Control"] <- "Infected"
df3$group[df3$group != "Control"] <- "Infected"
df4$group[df4$group != "Control"] <- "Infected"

df1$Study.Day <- df1$Day + hms_to_decimal_day(df1$chron_time)
df2$Study.Day <- df2$Day + hms_to_decimal_day(df2$chron_time)
df3$Study.Day <- df3$Day + hms_to_decimal_day(df3$chron_time)
df4$Study.Day <- df4$Day + hms_to_decimal_day(df4$chron_time)

df1$NHP <- "Cyno"
df1$virus <- "Dengue virus"
df2$NHP <- "Squirrel"
df2$virus <- "Dengue virus"
df3$NHP <- "Squirrel"
df3$virus <- "Zika virus"
df4$NHP <- "Cyno"
df4$virus <- "Zika virus"

# remove controls because already in another exp
df3 <- df3[df3$group != "Control",]
df4 <- df4[df4$group != "Control",]

df <- rbind(df1,df2,df3,df4)

df$monkey_status <- df$virus
df$monkey_status[which(df$group == "Control")] <- "Control"
df$group_agg <- interaction(df$NHP,df$monkey_status)
df$group_agg <- gsub(" ","_",df$group_agg) # no space for posthoc test to work
df$group_agg <- as.factor(df$group_agg)

# keep temperatures after infection has been performed
# and remove low temperatures, recorded after death
df <- df[df$Study.Day >= 0.5 & df$temp > 33,]

# Without day 28 
df <- df[df$Day != 28,]

m1 <- glmmTMB(temp ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # issues / expected with much data
testDispersion(m1) # ok
MuMIn::AICc(m1) # 298480.4

m2 <- glmmTMB(temp ~ group_agg + (1|ID) + (1|Day),
              data = df)
simulateResiduals(m2, plot = T) # issues / expected with much data
testDispersion(m2) # ok
MuMIn::AICc(m2) # 291424.7

anova(m1, m2) # in favour of random effects 

# select only the biologically relevant comparisons
tuk = glht(m2, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Zika_virus - Squirrel.Zika_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Cyno.Zika_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0")))
summary(tuk, test = adjusted("none")) # species diff all signif / idem
summary(tuk, test = adjusted("fdr")) # species diff all signif / idem

res <- tidy(confint(tuk))
# lower temp in cyno than squirrels





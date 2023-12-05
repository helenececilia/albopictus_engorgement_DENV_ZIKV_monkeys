## ---------------------------
##
## Script name:
##
## Purpose of script:
##
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

# Time metrics and weight (after day 0) / load data ----
df1 <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_with_total_correct.csv")
df1$NHP <- "Cyno"
df1$virus <- "Dengue virus"

df2 <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour_with_total_correct.csv")
df2$NHP <- "Squirrel"
df2$virus <- "Dengue virus"
df2$mosq_type <- NULL

df3 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour_with_total_correct.csv")
df3$NHP <- "Squirrel"
df3$virus <- "Zika virus"
df3$mosq_type <- NULL
# removing controls as they're already in df2
df3 <- df3[df3$group != "control",]

df1$raw_file_inspected <- NULL
df2$raw_file_inspected <- NULL
df3$raw_file_inspected <- NULL

df4 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Cyno_data_feeding_behaviour.csv")
df4$NHP <- "Cyno"
df4$virus <- "Zika virus"
df4$mosq_type <- NULL

df <- rbind(df1,df2,df3,df4)

df$monkey_status <- df$virus
df$monkey_status[which(df$group == "control")] <- "Control"

df$group_agg <- interaction(df$NHP,df$monkey_status)
df$group_agg <- gsub(" ","_",df$group_agg) # no space for posthoc test to work
df$group_agg <- as.factor(df$group_agg)

levels(df$group_agg)

feed_duration <- str_split(df$feed_duration,":")
vec_duration <- sapply(feed_duration,"[[",2)
df$int_duration <- as.integer(vec_duration)

df$is_day28 <- df$day == 28
df$is_day28 <- as.factor(df$is_day28)
df$time_of_day <- df$dec - df$day

# Test if day 28 different from other days ----
## Duration ----
summary(df$int_duration[df$is_day28 == FALSE])
summary(df$int_duration)

### Cyno ----
test <- df[df$NHP == "Cyno",]

summary(test$int_duration[test$is_day28 == TRUE])
summary(test$int_duration[test$is_day28 == FALSE])

m1 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = poisson) # for integer values, poisson accounting for overdispersion
# model convergence problem with nbinom1, so stick with poisson
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 562.6082 / 723.037 with ZIKV cyno

# m2 <- glmmTMB(int_duration ~ is_day28,
#               data = test,
#               family = nbinom1,
#               dispformula = ~ is_day28)
# model convergence problem

# m3 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
#               data = test,
#               family = nbinom1)
# model convergence problem

# m4 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
#               data = test,
#               dispformula = ~ is_day28,
#               family = nbinom1)
# model convergence problem

m2 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
              data = test,
              family = poisson) 
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # signif
MuMIn::AICc(m2) # 724.308 with ZIKV cyno

anova(m1,m2) # not signif

# m1 is the best we can do
plot(allEffects(m1, partial.residuals = T))
summary(m1) # no effect of day 28
car::Anova(m1, type = "II") # no effect of day 28

### Squirrel ----
test <- df[df$NHP == "Squirrel",]

summary(test$int_duration[test$is_day28 == TRUE])
summary(test$int_duration[test$is_day28 == FALSE])

# if you get this error : 
# Erreur dans exists(cacheKey, where = .rs.WorkingDataEnv, inherits = FALSE) : 
#   premier argument incorrect
# Erreur dans assign(cacheKey, frame, .rs.CachedDataEnv) : 
#   tentative d'utilisation de nom de variable de longueur nulle
# RESTART R (and reload your libraries)

m0 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = poisson)
simulateResiduals(m0, plot = T) # issues
testDispersion(m0) # signif
MuMIn::AICc(m0) # 698.1084

m1 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = nbinom1) # for integer values, poisson accounting for overdispersion
# false convergence
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 700.1834

anova(m0,m1) # not signif

m2 <- glmmTMB(int_duration ~ is_day28,
              data = test,
              family = nbinom1,
              dispformula = ~ is_day28)
# # model convergence problem

m3 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
              data = test,
              family = nbinom1)
# model convergence problem

m4 <- glmmTMB(int_duration ~ is_day28 + (1|ID),
              data = test,
              dispformula = ~ is_day28,
              family = nbinom1)
# model convergence problem

plot(allEffects(m0, partial.residuals = T))
summary(m0) # effect of day 28
car::Anova(m0, type = "II") # effect of day 28

# longer durations on day 28
res <- confint(m1)
exp(res) 

## Time of day ----
summary(df$time_of_day)
(0.6813*24*60)/60 # 16.3512
0.3512*60 # 21.1
summary(df$time_of_day[df$is_day28 == FALSE])
(0.3618056*24*60)/60 # 8.683334
0.683334*60 #41
(0.6038194*24*60)/60 # 14.49167
0.49167*60 # 29.5
df_no28 <- df[df$is_day28 == FALSE,]
hist(df$time_of_day)

### Cyno  ----
test <- df[df$NHP == "Cyno",]
m1 <- glmmTMB(time_of_day ~ is_day28,
              data = test) 
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
MuMIn::AICc(m1) # -564.0931 / -739.6126 with ZIKV cyno

m2 <- glmmTMB(time_of_day ~ is_day28 + (1|ID),
              data = test)
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
MuMIn::AICc(m2) # -662.3674 / -863.4266 with ZIKV cnyo

anova(m1,m2) # in favour of random effects

plot(allEffects(m2, partial.residuals = T))
summary(m2) # no effect of day 28
car::Anova(m2, type = "II") # no effect of day 28

### Squirrel ----
test <- df[df$NHP == "Squirrel",]
# if you get this error : 
# Erreur dans exists(cacheKey, where = .rs.WorkingDataEnv, inherits = FALSE) : 
#   premier argument incorrect
# Erreur dans assign(cacheKey, frame, .rs.CachedDataEnv) : 
#   tentative d'utilisation de nom de variable de longueur nulle
# RESTART R (and reload your libraries)

# Time of day is between 0 and 1, but 
# i) as it's a cycle, 0 and 1 are actually close in time 
# and ii) our values don't go near the extreme so we'll use linear models 
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

# Temperature : effect of species or infection -----

## All temperature data (every 15 mins) ----
df1 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Temperatures_Sylvatic_DENV-2_Cynomolgus_Macaques.csv")
df2 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Temperatures_Sylvatic_DENV-2_Squirrel_Monkeys.csv")
df3 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Temperatures_Sylvatic_ZIKV_Squirrel_Monkeys.csv")
df4 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Temperatures_Sylvatic_ZIKV_Cynomolgus_Macaques.csv")

treat1 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
                   dec = ".", sep  ="\t")
treat2 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat3 <- read.csv("~/Documents/POSTDOC/final_repositories/hanley_2023_sylvatic_DENV_ZIKV_trade_offs/data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv",
                   dec = ".", sep  ="\t")
treat4 <- read.csv("../data/ZIKV_Sylv_Cyno/master_ZIKV_Cyno.csv",
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

levels(df$group_agg)

# keep temperatures after infection has been performed
# and remove low temperatures, recorded after death
df <- df[df$Study.Day >= 0.5 & df$temp > 33,]
### With day 28 ----

# 2nd AICc values are with ZIKV cyno dataset
m1 <- glmmTMB(temp ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # issues / expected with much data?
testDispersion(m1) # ok
MuMIn::AICc(m1) # 285999.2 / 306139.8

m2 <- glmmTMB(temp ~ group_agg + (1|ID) + (1|Day),
              data = df)
simulateResiduals(m2, plot = T) # issues / expected with much data?
testDispersion(m2) # ok
MuMIn::AICc(m2) # 279429.3 / 298931.1

anova(m1, m2) # in favour of random effects / idem

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
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)
# lower temp in cyno than squirrels

### Without day 28 ----

df <- df[df$Day != 28,]

# 2nd AICc values are with ZIKV cyno dataset
m1 <- glmmTMB(temp ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # issues / expected with much data?
testDispersion(m1) # ok
MuMIn::AICc(m1) # 279011.2 / 298480.4

m2 <- glmmTMB(temp ~ group_agg + (1|ID) + (1|Day),
              data = df)
simulateResiduals(m2, plot = T) # issues / expected with much data?
testDispersion(m2) # ok
MuMIn::AICc(m2) # 272568.7 / 291424.7

anova(m1, m2) # in favour of random effects / idem

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
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)
# lower temp in cyno than squirrels

## Temperature at the time of feeding (after day 0) : effect of species or infection ----

df1 <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_with_total_correct.csv")
df1$NHP <- "Cyno"
df1$virus <- "Dengue virus"

df2 <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour_with_total_correct.csv")
df2$NHP <- "Squirrel"
df2$virus <- "Dengue virus"
df2$mosq_type <- NULL

df3 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour_with_total_correct.csv")
df3$NHP <- "Squirrel"
df3$virus <- "Zika virus"
df3$mosq_type <- NULL

df4 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Cyno_data_feeding_behaviour.csv")
df4$NHP <- "Cyno"
df4$virus <- "Zika virus"

df1$raw_file_inspected <- NULL
df2$raw_file_inspected <- NULL
df3$raw_file_inspected <- NULL

# removing controls if they're in other exp as well
df3 <- df3[df3$group != "control",]

df <- rbind(df1,df2,df3,df4)

df$monkey_status <- df$virus
df$monkey_status[which(df$group == "control")] <- "Control"

df$group_agg <- interaction(df$NHP,df$monkey_status)
df$group_agg <- gsub(" ","_",df$group_agg) # no space for posthoc test to work
df$group_agg <- as.factor(df$group_agg)

feed_duration <- str_split(df$feed_duration,":")
vec_duration <- sapply(feed_duration,"[[",2)
df$int_duration <- as.integer(vec_duration)

levels(df$group_agg)

### With day 28 ----

m1 <- glmmTMB(temp_estim_feed ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
MuMIn::AICc(m1) # 366.6067

m2 <- glmmTMB(temp_estim_feed ~ group_agg + (1|ID) + (1|day),
              data = df)
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
MuMIn::AICc(m2) # 321.2212

anova(m1, m2) # in favour of random effects

# select only the biologically relevant comparisons
tuk = glht(m2, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # species diff both signif + infected/control in cyno
summary(tuk, test = adjusted("fdr")) # species diff both signif + infected/control in cyno (close to 0.05)

plot(allEffects(m2, partial.residuals = T))
# higher temp in DENV-infected cynos than control, but CI contains 1 (odd ratio)
res <- tidy(confint(tuk))
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)

### Without day 28 ----

df <- df[df$day != 28,]

# 2nd AICc values are with ZIKV cyno
m1 <- glmmTMB(temp_estim_feed ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # ok
testDispersion(m1) # ok
MuMIn::AICc(m1) # 310.3752 / 350.8832

m2 <- glmmTMB(temp_estim_feed ~ group_agg + (1|ID) + (1|day),
              data = df)
simulateResiduals(m2, plot = T) # ok
testDispersion(m2) # ok
MuMIn::AICc(m2) # 266.8665 / 286.7474

anova(m1, m2) # in favour of random effects / idem

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
# higher temp in DENV-infected cynos than control, but CI contains 1 (odd ratio)
res <- tidy(confint(tuk))
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)


# Duration (after day 0) ----
## Effect of species and infection group / With day 28 ----
m1 <- glmmTMB(int_duration ~ group_agg,
              family = poisson,
              data = df)
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 1292.635

m2 <- glmmTMB(int_duration ~ group_agg,
              family = nbinom1,
              data = df)
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # signif
MuMIn::AICc(m2) # 1294.719

anova(m1,m2) # keep poisson

m3 <- glmmTMB(int_duration ~ group_agg + (1|ID) + (1|day),
              family = poisson,
              data = df)
simulateResiduals(m3, plot = T) # issues
testDispersion(m3) # signif
MuMIn::AICc(m3) # 1280.14

anova(m1, m3) # in favour of random effects

m4 <- glmmTMB(int_duration ~ group_agg + (1|ID) + (1|day),
              family = nbinom1,
              data = df)
# model convergence problem
simulateResiduals(m4, plot = T) # issues
testDispersion(m4) # signif
MuMIn::AICc(m4) # NA

anova(m3,m4) # not working

# select only the biologically relevant comparisons
tuk = glht(m3, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # effect of species, control only
summary(tuk, test = adjusted("fdr")) # nothing signif

plot(allEffects(m3, partial.residuals = T))
res <- tidy(confint(tuk))
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)
# lower duration in cyno control than squirrel control but CI contains 1

## Effect of species and infection group / Without day 28 ----
# 2nd AICc values are wtih ZIKV cyno
m1 <- glmmTMB(int_duration ~ group_agg,
              family = poisson,
              data = df[df$day != 28,])
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # signif
MuMIn::AICc(m1) # 1109.955 / 1250.184

m2 <- glmmTMB(int_duration ~ group_agg,
              family = nbinom1,
              data = df[df$day != 28,])
# model convergence problem
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # signif
MuMIn::AICc(m2) # 1112.051 / 1252.283

anova(m1,m2) # keep poisson / idem

m3 <- glmmTMB(int_duration ~ group_agg + (1|ID) + (1|day),
              family = poisson,
              data = df[df$day != 28,])
simulateResiduals(m3, plot = T) # issues
testDispersion(m3) # signif
MuMIn::AICc(m3) # 1112.163 / 1252.468

anova(m1, m3) # no need for random effects / idem

m4 <- glmmTMB(int_duration ~ group_agg + (1|ID) + (1|day),
              family = nbinom1,
              data = df[df$day != 28,])
# model convergence problem
simulateResiduals(m4, plot = T) # issues
testDispersion(m4) # signif
MuMIn::AICc(m4) # NA

anova(m3,m4) # not working

# select only the biologically relevant comparisons
tuk = glht(m1, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Zika_virus - Squirrel.Zika_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Cyno.Zika_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # nothing signif / longer durations for Cyno ZIKV than Cyno controls
summary(tuk, test = adjusted("fdr")) # nothing signif / longer durations for Cyno ZIKV than Cyno controls

plot(allEffects(m1, partial.residuals = T))
res <- tidy(confint(tuk))
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)

new_dat <- data.frame(group_agg = c("Cyno.Dengue_virus","Cyno.Control","Squirrel.Dengue_virus",
                            "Squirrel.Control","Squirrel.Zika_virus","Cyno.Zika_virus"))

pp <- predict(m1, se.fit = TRUE, newdata = new_dat)
# the inverse link of a poisson is exp
ci_lwr <- with(pp, exp(fit + qnorm(0.025)*se.fit))
ci_upr <- with(pp, exp(fit + qnorm(0.975)*se.fit))
means <- with(pp, exp(fit))
names(means) <- new_dat$group
names(ci_upr) <- new_dat$group
names(ci_lwr) <- new_dat$group
means
ci_lwr
ci_upr

# Weight (after day 0) ----
## Effect of species and infection group / With day 28 ----
m1 <- glmmTMB(weight ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # 4187.286

m2 <- glmmTMB(weight ~ group_agg + (1|ID) + (1|day),
              data = df)
# model convergence problem
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # ok
MuMIn::AICc(m2) # NA

anova(m1, m2) # not working

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

## Effect of species and infection group / Without day 28 ----

m1 <- glmmTMB(weight ~ group_agg,
              data = df[df$day != 28,])
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # 3701.645

m2 <- glmmTMB(weight ~ group_agg + (1|ID) + (1|day),
              data = df[df$day != 28,])
# error in eigen(h)

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

# Time of day (after day 0) : DO NOT INCLUDE / NOT RELEVANT? ----

## Effect of species and infection group / With day 28 ----
levels(df$group_agg)

m1 <- glmmTMB(time_of_day ~ group_agg,
              data = df)
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # -878.9825

m2 <- glmmTMB(time_of_day ~ group_agg + (1|ID) + (1|day),
              data = df)
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # ok
MuMIn::AICc(m2) # -927.3261

anova(m1, m2) # in favour of random effects

# select only the biologically relevant comparisons
tuk = glht(m2, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # effect of infection (all)
summary(tuk, test = adjusted("fdr")) # effect of infection (all)
# due to sequential mosquito feeding / controls before infected

plot(allEffects(m2, partial.residuals = T))
res <- tidy(confint(tuk))

## Effect of species and infection group / Without day 28 ----

m1 <- glmmTMB(time_of_day ~ group_agg,
              data = df[df$day != 28,])
simulateResiduals(m1, plot = T) # issues
testDispersion(m1) # ok
MuMIn::AICc(m1) # -967.5297

m2 <- glmmTMB(time_of_day ~ group_agg + (1|ID) + (1|day),
              data = df[df$day != 28,])
simulateResiduals(m2, plot = T) # issues
testDispersion(m2) # ok
MuMIn::AICc(m2) # -1007.434

anova(m1, m2) # in favour of random effects

# select only the biologically relevant comparisons
tuk = glht(m2, linfct = mcp(group_agg = c("Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
                                          "Cyno.Control - Squirrel.Control = 0",
                                          "Cyno.Dengue_virus - Cyno.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Control = 0",
                                          "Squirrel.Zika_virus - Squirrel.Control = 0",
                                          "Squirrel.Dengue_virus - Squirrel.Zika_virus = 0"))) 
summary(tuk, test = adjusted("none")) # effect of infection (ZIKV squirrel and DENv cyno) + species diff DENV infected + DENV/ZIKV squirrel diff
summary(tuk, test = adjusted("fdr")) # idem
# due to sequential mosquito feeding / controls before infected

plot(allEffects(m2, partial.residuals = T))
res <- tidy(confint(tuk))

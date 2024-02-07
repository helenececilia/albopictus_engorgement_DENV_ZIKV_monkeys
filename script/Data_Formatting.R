## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Helene Cecilia
##
## Date Created: 2023-02-20

rm(list=ls())

## Loading Packages  ------------------
require(ggdist)
require(tidyverse)
require(patchwork)
require(Matching) # for ks.boot
require(scales)
library(DHARMa)
library(glmmTMB)
library(effects)
library(MASS)
library(mgcv)
library(readxl)
library(lubridate)
library(hms)
library(chron)
library(MuMIn) # for r.squaredGLMM
# library(gammit) # for predict_gamm
library(janitor) # for clean_names

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # print data frame
  return(data_frame)
}


hms_to_decimal_day <-function(time_hms){
  h <- hours(time_hms)
  m <- minutes(time_hms)
  s <- seconds(time_hms)
  sum_minutes <- h*60 + m + s/60
  dec <- sum_minutes / (24*60)
  return(dec)
}

decimal_day_to_hms <- function(time_dec){
  float_hours <- time_dec * 24
  int_hours <- floor(float_hours)
  int_minutes <- as.integer((float_hours - int_hours)*60) # we don't bother with seconds
  if(int_minutes < 10){
    int_minutes <- paste0("0",int_minutes)
  }
  str_time <- paste(int_hours,int_minutes,"00",sep = ":")
  return(str_time)
  # time_hms <- as_hms(str_time)
  # browser()
  # return(time_hms)
}

# Returns the indices of minimum values
which.mins <- function(x, mins=2) {
  head(order(x), mins)
}

# ZIKV cyno TEST ----

# New temperature data ----
# Feeding time data ----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_zikv_cy <- feeding_time$`ZIKV Cyno`
feeding_zikv_cy <- feeding_zikv_cy[,seq(1,6)]
  
test <- feeding_zikv_cy %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                   date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

# # descriptive stats - duration
# as_hms(mean(test$diff[test$Day == 0]))
# as_hms(median(test$diff[test$Day == 0]))
# as_hms(min(test$diff[test$Day == 0]))
# as_hms(max(test$diff[test$Day == 0]))
# 
# # with/without outliers
# as_hms(mean(test$diff[test$Day != 0]))
# as_hms(median(test$diff[test$Day != 0]))
# as_hms(min(test$diff[test$Day != 0]))
# as_hms(max(test$diff[test$Day != 0]))
# 
# as_hms(mean(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(median(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(min(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(max(test$diff[test$Day != 0 & is.na(test$outlier)]))
# 
# # descriptive stats - time of day
# as_hms(mean(test$mid_time[test$Day == 0]))
# as_hms(median(test$mid_time[test$Day == 0]))
# as_hms(min(test$mid_time[test$Day == 0]))
# as_hms(max(test$mid_time[test$Day == 0]))
# 
# # with/without outliers
# as_hms(mean(test$mid_time[test$Day != 0]))
# as_hms(median(test$mid_time[test$Day != 0]))
# as_hms(min(test$mid_time[test$Day != 0]))
# as_hms(max(test$mid_time[test$Day != 0]))
# 
# as_hms(mean(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(median(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(min(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(max(test$mid_time[test$Day != 0 & is.na(test$outlier)]))


test$diff <- as_hms(test$diff)
feed_zikv_cy_select <- test[,c("NHP","Day","mid_time","diff","outlier")]
colnames(feed_zikv_cy_select) <- c("ID","Day","feed_mean","feed_duration","outlier")

# Spaces in ID ; matching will not work
feed_zikv_cy_select$ID <- gsub(" ","",feed_zikv_cy_select$ID, fixed = TRUE)

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ID", "Day", "temp_estim_feed"))
for(i in unique(feed_zikv_cy_select$ID)){
  temp_df <- read.table(paste0("../data/monkey_temp_new_names/Cyno_",i,".txt"),
                        dec = ".", sep = "\t")
  temp_df <- temp_df[-1,] # remove first line (headers that could not be read correctly)
  colnames(temp_df) <- c("X","temp")
  temp_df$date <- mdy_hm(temp_df$X) # adapt format here if needed
  # midnight times are not explicit (only mdy in X)
  # returns NA after conversion, correct manually
  na_date <- temp_df$X[which(is.na(temp_df$date))]
  na_date <- paste(na_date,"0:00", sep = " ")
  temp_df$date[which(is.na(temp_df$date))] <- mdy_hm(na_date) # adapt format here if needed
  
  # compute days since/before infection
  temp_df$Day <- difftime(as.Date(temp_df$date), as.Date("2022-05-02"), units = "days")

  temp_df$Day <- as.integer(temp_df$Day)
  temp_df$X <- NULL
  temp_df$temp <- as.numeric(temp_df$temp)
  
  for(d in unique(feed_zikv_cy_select$Day)){
    day_monk <- feed_zikv_cy_select[feed_zikv_cy_select$ID == i & feed_zikv_cy_select$Day == d,]
    feed_mean <- unique(day_monk$feed_mean)
    if(dim(day_monk)[[1]] != 0){ #!is.na(feed_mean)
      temp_df_day <- temp_df[temp_df$Day == d,]
      temp_df_day$time <- as_hms(temp_df_day$date)
      
      df <- merge(day_monk,temp_df_day, by = c("Day"))
      # browser()
      # find in which window to do linear interpolation
      select <- df[which.mins(abs(df$time - df$feed_mean)),]
      select <- select[order(select$time),] # increasing times
      # browser()
      temp_diff <- select$temp[2] - select$temp[1]
      time_diff <- as.numeric(select$time[2] - select$time[1])
      slope <- temp_diff/time_diff
      temp_estim <- select$temp[1] +
        slope * as.numeric(feed_mean - select$time[1])
      res <- rbind(res,data.frame(ID = i,
                                  Day = d,
                                  temp_estim_feed = temp_estim))
      # browser()
    }
  }
}
new <- merge(feed_zikv_cy_select, res, by = c("ID","Day"))
end <- Sys.time()
end-start # 0.2 scd

write.csv(new, "../output/mosq_feeding_behaviour/ZIKV_Cyno_temperature_estimate_at_feeding.csv", row.names = F)

## Converting feeding times into day and its decimal ----
zikv_cy_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Cyno_temperature_estimate_at_feeding.csv")

zikv_cy_feed_time_and_temp$chron_time <- chron(times. = as.character(zikv_cy_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

zikv_cy_feed_time_and_temp$dec <- zikv_cy_feed_time_and_temp$Day + hms_to_decimal_day(zikv_cy_feed_time_and_temp$chron_time)

## Including feeding behaviour ----
colnames(zikv_cy_feed_time_and_temp)[2] <- "day"

select_temp_cyno <- unique(zikv_cy_feed_time_and_temp[,c("ID","day","temp_estim_feed","dec","feed_duration","outlier")])

df <- read.csv("../data/ZIKV_Sylv_Cyno/master_ZIKV_Cyno.csv",
               sep = "\t", dec = ".") 
# TEMPORARY : NO DEDUCED YET, USE RAW VIREMIA
df$Viremia.deduced <- df$Viremia..log10.pfu.ml..1.

if(max(df$Viremia.deduced, na.rm = T) >= 10){ # then viremia are not in log10 yet
  df$Viremia.deduced <- log10(df$Viremia.deduced+1) 
}
df$group <- df$Final.Treatment
df <- df[,c("ID","Day.Post.Infection","Viremia.deduced",
            "Number.of.IT.injected.mosquitoes.that.fed",
            "Number.of.mosquitoes.that.engorged","group","Sex","Weight..g.",
            "Total.number.of.mosquitoes.exposed..sometimes.differs.from.the.expected.15.")]
colnames(df) <- c("ID","day","viremia",
                  "nb_fed_day0","nb_fed",
                  "group","sex","weight",
                  "total_correct")
if(max(df$weight, na.rm = T) < 50){ # then weights are not in grams
  df$weight <- 1000*df$weight # to be in grams
}
df$group[df$group != "Control"] <- "infected"
df$group <- tolower(df$group) # remove upper case for Control

# Dataset day 0
day_0 <- df[df$day == 0,]
day_0 <- day_0[,c("ID","day","nb_fed_day0","group","weight","sex","total_correct")]
day_0$total <- day_0$total_correct #ifelse(day_0$nb_fed_day0 > 1, 10, 1)
day_0_with_temp <- merge(day_0, select_temp_cyno, by = c("ID","day"))
write.csv(day_0_with_temp, "../output/mosq_feeding_behaviour/ZIKV_Cyno_data_feeding_behaviour_day0.csv",
          row.names = F)

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
# others <- others[complete.cases(others),] # PUT AGAIN ONCE WE HAVE SEX
zikv_cy <- merge(others, select_temp_cyno, by = c("ID","day"))

zikv_cy$total <- zikv_cy$total_correct #ifelse(zikv_cy$nb_fed <= 10, 10, zikv_cy$nb_fed)

write.csv(zikv_cy, "../output/mosq_feeding_behaviour/ZIKV_Cyno_data_feeding_behaviour.csv",
          row.names = F)


# DENV cyno ----

# New temperature data ----
# Feeding time data ----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_denv_cy <- feeding_time$`DENV Cyno`
feeding_denv_cy <- feeding_denv_cy[,seq(1,6)]
test <- feeding_denv_cy %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                   date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

# # descriptive stats - duration
# as_hms(mean(test$diff[test$Day == 0]))
# as_hms(median(test$diff[test$Day == 0]))
# as_hms(min(test$diff[test$Day == 0]))
# as_hms(max(test$diff[test$Day == 0]))
# 
# # with/without outliers
# as_hms(mean(test$diff[test$Day != 0]))
# as_hms(median(test$diff[test$Day != 0]))
# as_hms(min(test$diff[test$Day != 0]))
# as_hms(max(test$diff[test$Day != 0]))
# 
# as_hms(mean(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(median(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(min(test$diff[test$Day != 0 & is.na(test$outlier)]))
# as_hms(max(test$diff[test$Day != 0 & is.na(test$outlier)]))
# 
# # descriptive stats - time of day
# as_hms(mean(test$mid_time[test$Day == 0]))
# as_hms(median(test$mid_time[test$Day == 0]))
# as_hms(min(test$mid_time[test$Day == 0]))
# as_hms(max(test$mid_time[test$Day == 0]))
# 
# # with/without outliers
# as_hms(mean(test$mid_time[test$Day != 0]))
# as_hms(median(test$mid_time[test$Day != 0]))
# as_hms(min(test$mid_time[test$Day != 0]))
# as_hms(max(test$mid_time[test$Day != 0]))
# 
# as_hms(mean(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(median(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(min(test$mid_time[test$Day != 0 & is.na(test$outlier)]))
# as_hms(max(test$mid_time[test$Day != 0 & is.na(test$outlier)]))


test$diff <- as_hms(test$diff)
feed_denv_cy_select <- test[,c("NHP","Day","mid_time","diff","outlier")]
colnames(feed_denv_cy_select) <- c("ID","Day","feed_mean","feed_duration","outlier")

# Spaces in ID ; matching will not work
feed_denv_cy_select$ID <- gsub(" ","",feed_denv_cy_select$ID, fixed = TRUE)

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ID", "Day", "temp_estim_feed"))
for(i in unique(feed_denv_cy_select$ID)){
  
  temp_df <- read.table(paste0("../data/monkey_temp_new_names/Cyno_",i,".txt"),
                        dec = ".", sep = "\t")
  temp_df <- temp_df[-1,] # remove first line (headers that could not be read correctly)
  colnames(temp_df) <- c("X","temp")
  temp_df$date <- mdy_hms(temp_df$X)
  # midnight times are not explicit (only mdy in X)
  # returns NA after conversion, correct manually
  na_date <- temp_df$X[which(is.na(temp_df$date))]
  na_date <- paste(na_date,"12:00:00 AM", sep = " ")
  temp_df$date[which(is.na(temp_df$date))] <- mdy_hms(na_date)
  if(i %notin% c("UG253A","BC407")){
    temp_df$Day <- difftime(as.Date(temp_df$date), as.Date("2020-11-09"), units = "days")
  }else{
    temp_df$Day <- difftime(as.Date(temp_df$date), as.Date("2020-11-10"), units = "days")
  }
  temp_df$Day <- as.integer(temp_df$Day)
  temp_df$X <- NULL
  temp_df$temp <- as.numeric(temp_df$temp)
  
  for(d in unique(feed_denv_cy_select$Day)){
    day_monk <- feed_denv_cy_select[feed_denv_cy_select$ID == i & feed_denv_cy_select$Day == d,]
    feed_mean <- unique(day_monk$feed_mean)
    if(!is.na(feed_mean)){
      temp_df_day <- temp_df[temp_df$Day == d,]
      temp_df_day$time <- as_hms(temp_df_day$date)
      
      df <- merge(day_monk,temp_df_day, by = c("Day"))
      # browser()
      # find in which window to do linear interpolation
      select <- df[which.mins(abs(df$time - df$feed_mean)),]
      select <- select[order(select$time),] # increasing times
      # browser()
      temp_diff <- select$temp[2] - select$temp[1]
      time_diff <- as.numeric(select$time[2] - select$time[1])
      slope <- temp_diff/time_diff
      temp_estim <- select$temp[1] +
        slope * as.numeric(feed_mean - select$time[1])
      res <- rbind(res,data.frame(ID = i,
                                  Day = d,
                                  temp_estim_feed = temp_estim))
      # browser()
    }
  }
}
new <- merge(feed_denv_cy_select, res, by = c("ID","Day"))
end <- Sys.time()
end-start # 1.8 scd

write.csv(new, "../output/mosq_feeding_behaviour/DENV_Cyno_temperature_estimate_at_feeding.csv", row.names = F)

## Converting feeding times into day and its decimal ----
denv_cy_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_temperature_estimate_at_feeding.csv")

denv_cy_feed_time_and_temp$chron_time <- chron(times. = as.character(denv_cy_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

denv_cy_feed_time_and_temp$dec <- denv_cy_feed_time_and_temp$Day + hms_to_decimal_day(denv_cy_feed_time_and_temp$chron_time)

## Including feeding behaviour ----
colnames(denv_cy_feed_time_and_temp)[2] <- "day"

select_temp_cyno <- unique(denv_cy_feed_time_and_temp[,c("ID","day","temp_estim_feed","dec","feed_duration","outlier")])

df <- read.csv("../data/Table_S1_Sylvatic_DENV-2_Cynomolgus_Macaques.csv",
               sep = "\t", dec = ".")

if(max(df$Viremia.deduced, na.rm = T) >= 10){ # then viremia are not in log10 yet
  df$Viremia.deduced <- log10(df$Viremia.deduced+1) 
}
df$group <- df$Final.Treatment
df <- df[,c("ID","Day.Post.Infection","Viremia.deduced",
            "Number.of.IT.injected.mosquitoes.that.fed",
            "Number.of.mosquitoes.that.engorged","group","Sex","Weight..kg.")]
colnames(df) <- c("ID","day","viremia",
                  "nb_fed_day0","nb_fed",
                  "group","sex","weight")
if(max(df$weight, na.rm = T) < 50){ # then weights are not in grams
  df$weight <- 1000*df$weight # to be in grams
}
df$group[df$group != "Control"] <- "infected"
df$group <- tolower(df$group) # remove upper case for Control

# Dataset day 0
day_0 <- df[df$day == 0,]
day_0 <- day_0[,c("ID","day","nb_fed_day0","group","weight","sex")]
day_0_with_temp <- merge(day_0, select_temp_cyno, by = c("ID","day"))
day_0_with_temp$total <- ifelse(day_0_with_temp$nb_fed_day0 > 1, 10, 1)

write.csv(day_0_with_temp, "../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_day0.csv",
          row.names = F)

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
others <- others[complete.cases(others),] 
denv_cy <- merge(others, select_temp_cyno, by = c("ID","day"))

# 10 is the default, occasionnally more
denv_cy$total <- ifelse(denv_cy$nb_fed <= 10, 10, denv_cy$nb_fed)

write.csv(denv_cy, "../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour.csv",
          row.names = F)


# ZIKV squirrel ------
# New temperature data ----
# Feeding time data ----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_zikv <- feeding_time$`ZIKV Saimiri`
feeding_zikv$Notes <- NULL
test <- feeding_zikv %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

# # descriptive stats - duration
# as_hms(mean(test$diff[test$Day == 0],na.rm = T))
# as_hms(median(test$diff[test$Day == 0],na.rm = T))
# as_hms(min(test$diff[test$Day == 0],na.rm = T))
# as_hms(max(test$diff[test$Day == 0],na.rm = T))
# # with/without outliers
# as_hms(mean(test$diff[test$Day != 0],na.rm = T))
# as_hms(median(test$diff[test$Day != 0],na.rm = T))
# as_hms(min(test$diff[test$Day != 0],na.rm = T))
# as_hms(max(test$diff[test$Day != 0],na.rm = T))
# 
# as_hms(mean(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(median(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(min(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(max(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# 
# # descriptive stats - time of day
# as_hms(mean(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(median(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(min(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(max(test$mid_time[test$Day == 0],na.rm = T))
# 
# # with/without outliers
# as_hms(mean(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(median(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(min(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(max(test$mid_time[test$Day != 0],na.rm = T))
# 
# as_hms(mean(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(median(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(min(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(max(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))

test$diff <- as_hms(test$diff)
feed_zikv_select <- test[,c("NHP","Day","mid_time","mosq_type","diff","outlier")]
colnames(feed_zikv_select) <- c("ID","Day","feed_mean","mosq_type","feed_duration","outlier")
feed_zikv_select <- feed_zikv_select[complete.cases(feed_zikv_select[,c("feed_mean",
                                                                        "mosq_type",
                                                                        "feed_duration")]),]

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("ID", "Day",
                                                          "temp_estim_feed",
                                                          "mosq_type"))
for(i in unique(feed_zikv_select$ID)){
  
  temp_df <- read.table(paste0("../data/monkey_temp_new_names/Saimiri_",i,".txt"),
                        dec = ".", sep = "\t")
  temp_df <- temp_df[-1,] # remove first line (headers that could not be read correctly)
  colnames(temp_df) <- c("X","temp")
  temp_df$date <- mdy_hms(temp_df$X)
  # midnight times are not explicit (only mdy in X)
  # returns NA after conversion, correct manually
  na_date <- temp_df$X[which(is.na(temp_df$date))]
  na_date <- paste(na_date,"12:00:00 AM", sep = " ")
  temp_df$date[which(is.na(temp_df$date))] <- mdy_hms(na_date)

  temp_df$Day <- difftime(as.Date(temp_df$date), as.Date("2021-04-12"), units = "days")

  temp_df$Day <- as.integer(temp_df$Day)
  temp_df$X <- NULL
  temp_df$temp <- as.numeric(temp_df$temp)
  
  for(d in unique(feed_zikv_select$Day)){
    day_monk <- feed_zikv_select[feed_zikv_select$ID == i & feed_zikv_select$Day == d,]
    temp_df_day <- temp_df[temp_df$Day == d,]
    temp_df_day$time <- as_hms(temp_df_day$date)
    for(t in unique(day_monk$mosq_type)){ # for day 0, 2 different types, each with their own computation
      test_t <- day_monk[day_monk$mosq_type == t,]
      feed_mean <- unique(test_t$feed_mean)
      if(!is.na(feed_mean)){
        df <- merge(test_t,temp_df_day, by = c("Day"))
        # browser()
        # find in which window to do linear interpolation
        select <- df[which.mins(abs(df$time - df$feed_mean)),]
        select <- select[order(select$time),] # increasing times
        # browser()
        temp_diff <- select$temp[2] - select$temp[1]
        time_diff <- as.numeric(select$time[2] - select$time[1])
        slope <- temp_diff/time_diff
        temp_estim <- select$temp[1] +
          slope * as.numeric(feed_mean - select$time[1])
        res <- rbind(res,data.frame(ID = i,
                                    Day = d,
                                    temp_estim_feed = temp_estim,
                                    mosq_type = t))
        # browser()
      }
    }
  }
}
new <- merge(feed_zikv_select, res, by = c("ID","Day","mosq_type"))
end <- Sys.time()
end-start # 0.86 scd

write.csv(new, "../output/mosq_feeding_behaviour/ZIKV_Squirrel_temperature_estimate_at_feeding.csv",
          row.names = F)


## Converting feeding times into day and its decimal ----
zikv_sq_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_temperature_estimate_at_feeding.csv")

zikv_sq_feed_time_and_temp$chron_time <- chron(times. = as.character(zikv_sq_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

zikv_sq_feed_time_and_temp$dec <- zikv_sq_feed_time_and_temp$Day + hms_to_decimal_day(zikv_sq_feed_time_and_temp$chron_time)


# DENV squirrel -----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_denv_sq <- feeding_time$`DENV Saimiri`
feeding_denv_sq$Notes <- NULL
test <- feeding_denv_sq %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                   date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

# # descriptive stats - duration
# as_hms(mean(test$diff[test$Day == 0],na.rm = T))
# as_hms(median(test$diff[test$Day == 0],na.rm = T))
# as_hms(min(test$diff[test$Day == 0],na.rm = T))
# as_hms(max(test$diff[test$Day == 0],na.rm = T))
# # with/without outliers
# as_hms(mean(test$diff[test$Day != 0],na.rm = T))
# as_hms(median(test$diff[test$Day != 0],na.rm = T))
# as_hms(min(test$diff[test$Day != 0],na.rm = T))
# as_hms(max(test$diff[test$Day != 0],na.rm = T))
# 
# as_hms(mean(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(median(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(min(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(max(test$diff[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# 
# # descriptive stats - time of day
# as_hms(mean(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(median(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(min(test$mid_time[test$Day == 0],na.rm = T))
# as_hms(max(test$mid_time[test$Day == 0],na.rm = T))
# 
# # with/without outliers
# as_hms(mean(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(median(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(min(test$mid_time[test$Day != 0],na.rm = T))
# as_hms(max(test$mid_time[test$Day != 0],na.rm = T))
# 
# as_hms(mean(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(median(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(min(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))
# as_hms(max(test$mid_time[test$Day != 0 & is.na(test$outlier)],na.rm = T))

test$diff <- as_hms(test$diff)
feed_denv_sq_select <- test[,c("NHP","Day","mid_time","mosq_type","diff","outlier")]
colnames(feed_denv_sq_select) <- c("ID","Day","feed_mean","mosq_type","feed_duration","outlier")
feed_denv_sq_select <- feed_denv_sq_select[complete.cases(feed_denv_sq_select[,c("feed_mean",
                                                                                 "mosq_type",
                                                                                 "feed_duration")]),]

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("ID", "Day",
                                                          "temp_estim_feed",
                                                          "mosq_type"))
for(i in unique(feed_denv_sq_select$ID)){
  
  temp_df <- read.table(paste0("../data/monkey_temp_new_names/Saimiri_",i,".txt"),
                        dec = ".", sep = "\t")
  temp_df <- temp_df[-1,] # remove first line (headers that could not be read correctly)
  colnames(temp_df) <- c("X","temp")
  temp_df$date <- mdy_hms(temp_df$X)
  # midnight times are not explicit (only mdy in X)
  # returns NA after conversion, correct manually
  na_date <- temp_df$X[which(is.na(temp_df$date))]
  na_date <- paste(na_date,"12:00:00 AM", sep = " ")
  temp_df$date[which(is.na(temp_df$date))] <- mdy_hms(na_date)
  
  temp_df$Day <- difftime(as.Date(temp_df$date), as.Date("2021-04-05"), units = "days")
  
  temp_df$Day <- as.integer(temp_df$Day)
  temp_df$X <- NULL
  temp_df$temp <- as.numeric(temp_df$temp)
  
  for(d in unique(feed_denv_sq_select$Day)){
    day_monk <- feed_denv_sq_select[feed_denv_sq_select$ID == i & feed_denv_sq_select$Day == d,]
    temp_df_day <- temp_df[temp_df$Day == d,]
    temp_df_day$time <- as_hms(temp_df_day$date)
    for(t in unique(day_monk$mosq_type)){ # for day 0, 2 different types, each with their own computation
      test_t <- day_monk[day_monk$mosq_type == t,]
      feed_mean <- unique(test_t$feed_mean)
      if(!is.na(feed_mean)){
        df <- merge(test_t,temp_df_day, by = c("Day"))
        # browser()
        # find in which window to do linear interpolation
        select <- df[which.mins(abs(df$time - df$feed_mean)),]
        select <- select[order(select$time),] # increasing times
        # browser()
        temp_diff <- select$temp[2] - select$temp[1]
        time_diff <- as.numeric(select$time[2] - select$time[1])
        slope <- temp_diff/time_diff
        temp_estim <- select$temp[1] +
          slope * as.numeric(feed_mean - select$time[1])
        res <- rbind(res,data.frame(ID = i,
                                    Day = d,
                                    temp_estim_feed = temp_estim,
                                    mosq_type = t))
        # browser()
      }
    }
  }
}
new <- merge(feed_denv_sq_select, res, by = c("ID","Day","mosq_type"))
end <- Sys.time()
end-start # 0.9 scd

write.csv(new, "../output/mosq_feeding_behaviour/DENV_Squirrel_temperature_estimate_at_feeding.csv",
          row.names = F)

## Converting feeding times into day and its decimal ----
denv_sq_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_temperature_estimate_at_feeding.csv")

denv_sq_feed_time_and_temp$chron_time <- chron(times. = as.character(denv_sq_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

denv_sq_feed_time_and_temp$dec <- denv_sq_feed_time_and_temp$Day + hms_to_decimal_day(denv_sq_feed_time_and_temp$chron_time)

# Including feeding behaviour ----
# Aggregated controls squirrels in zikv and denv datasets
colnames(zikv_sq_feed_time_and_temp)[2] <- "day"
colnames(denv_sq_feed_time_and_temp)[2] <- "day"

temp_squirrel <- rbind(zikv_sq_feed_time_and_temp, denv_sq_feed_time_and_temp)
select_temp_squirrel <- unique(temp_squirrel[,c("ID","day","mosq_type","temp_estim_feed","dec",
                                                "feed_duration","outlier")])

df <- read.csv("../data/ZIKV_Sylv_Squirrel/Trade-Off_Project_Saimiri_Inf_Sylvatic_ZIKV_2021_Comprehensive_data_14Mar23_HC.csv", dec = ".", sep = "\t")
df <- df[df$Final.Treatment != "Control",]
if(max(df$Viremia.deduced, na.rm = T) >= 10){ # then viremia are not in log10 yet
  df$Viremia.deduced <- log10(df$Viremia.deduced+1) 
}
df$group <- df$Final.Treatment
df <- df[,c("ID","Day.Post.Infection","Viremia.deduced",
            "Number.of.IT.injected.mosquitoes.that.fed",
            "Number.of.mosquitoes.that.engorged","group","Sex","Weight..g.")]
colnames(df) <- c("ID","day","viremia",
                  "nb_fed_day0","nb_fed",
                  "group","sex","weight")
if(max(df$weight, na.rm = T) < 50){ # then weights are not in grams
  df$weight <- 1000*df$weight # to be in grams
}
df$group <- "infected"

# Dataset day 0
day_0 <- df[df$day == 0,]
day_0 <- day_0[,c("ID","day","nb_fed_day0","group","weight","sex")]
day_0$total <- 15
day_0_inf <- day_0
day_0_inf$mosq_type <- "infecting"
day_0_cof <- day_0
day_0_cof$mosq_type <- "cofeeding"
day_0_cof$nb_fed_day0 <- c(8,15,9,10,3,11,6,14,4,12)
day_0 <- rbind(day_0_inf,day_0_cof)
day_0_with_temp_zikv_inf <- merge(day_0, select_temp_squirrel, by = c("ID","day","mosq_type"))

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
others <- others[complete.cases(others),] 
zikv_inf_others <- merge(others, select_temp_squirrel, by = c("ID","day"))


df <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv",
               dec = ".", sep = "\t")
df <- df[df$Final.Treatment != "Control",]
if(max(df$Viremia.deduced, na.rm = T) >= 10){ # then viremia are not in log10 yet
  df$Viremia.deduced <- log10(df$Viremia.deduced+1) 
}
df$group <- df$Final.Treatment
df <- df[,c("ID","Day.Post.Infection","Viremia.deduced",
            "Number.of.IT.injected.mosquitoes.that.fed",
            "Number.of.mosquitoes.that.engorged","group","Sex","Weight..g.")]
colnames(df) <- c("ID","day","viremia",
                  "nb_fed_day0","nb_fed",
                  "group","sex","weight")
if(max(df$weight, na.rm = T) < 50){ # then weights are not in grams
  df$weight <- 1000*df$weight # to be in grams
}
df$group <- "infected"

# Dataset day 0
day_0 <- df[df$day == 0,]
day_0 <- day_0[,c("ID","day","nb_fed_day0","group","weight","sex")]
day_0$total <- 15
day_0_inf <- day_0
day_0_inf$mosq_type <- "infecting"
day_0_cof <- day_0
day_0_cof$mosq_type <- "cofeeding"
day_0_cof$nb_fed_day0 <- c(4,7,8,12,13,8,10,2,7,10)
day_0 <- rbind(day_0_inf,day_0_cof)
day_0_with_temp_denv_inf <- merge(day_0, select_temp_squirrel, by = c("ID","day","mosq_type"))

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
others <- others[complete.cases(others),] 
denv_inf_others <- merge(others, select_temp_squirrel, by = c("ID","day"))

controls_zikv_sq <- read.csv("../data/Table_S3_Sylvatic_ZIKV_Squirrel_Monkeys.csv", dec = ".", sep = "\t")
controls_zikv_sq <- controls_zikv_sq[controls_zikv_sq$Final.Treatment == "Control",]

controls_denv_sq <- read.csv("../data/Table_S2_Sylvatic_DENV-2_Squirrel_Monkeys.csv", dec = ".", sep = "\t")
controls_denv_sq <- controls_denv_sq[controls_denv_sq$Final.Treatment == "Control",]

controls <- rbind(controls_denv_sq[,c("ID","Day.Post.Infection",
                                      "Number.of.IT.injected.mosquitoes.that.fed",
                                      "Number.of.mosquitoes.that.engorged",
                                      "Sex","Weight..g.")],
                  controls_zikv_sq[,c("ID","Day.Post.Infection",
                                      "Number.of.IT.injected.mosquitoes.that.fed",
                                      "Number.of.mosquitoes.that.engorged",
                                      "Sex","Weight..g.")])
colnames(controls) <- c("ID","day","nb_fed_day0","nb_fed","sex","weight")
if(max(controls$weight, na.rm = T) < 50){ # then weights are not in grams
  controls$weight <- 1000*controls$weight # to be in grams
}
controls$group <- "control"

day_0 <- controls[controls$day == 0,]
day_0 <- day_0[,c("ID","day","nb_fed_day0","group","weight","sex")]
day_0$total <- 15

controls_day0 <- merge(day_0, select_temp_squirrel, by = c("ID","day"))
zikv_day0 <- rbind(day_0_with_temp_zikv_inf, controls_day0)
denv_day0 <- rbind(day_0_with_temp_denv_inf, controls_day0)

write.csv(zikv_day0, "../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour_day0.csv",
          row.names = F)
write.csv(denv_day0, "../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour_day0.csv",
          row.names = F)

control_others <- controls[controls$day > 0,]
control_others$nb_fed_day0 <- NULL
control_others$viremia <- 0
control_others <- control_others[complete.cases(control_others),]
control_others <- merge(control_others, select_temp_squirrel, by = c("ID","day"))

zikv_data <- rbind(zikv_inf_others, control_others)
zikv_data$total <- 15

denv_data <- rbind(denv_inf_others, control_others)
denv_data$total <- 15

write.csv(zikv_data, "../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour.csv", row.names = F)
write.csv(denv_data, "../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour.csv", row.names = F)

# Descriptive stats global dataset -----

# df1 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour_day0.csv")
# df2 <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour_day0.csv")
# df3 <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_day0.csv")
# df1$mosq_type <- NULL
# df2$mosq_type <- NULL
# df <- rbind(df1,df2,df3)
# df <- unique(df) # remove duplicate controls

df1 <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_data_feeding_behaviour.csv")
df2 <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_data_feeding_behaviour.csv")
df3 <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour.csv")
df1$mosq_type <- NULL
df2$mosq_type <- NULL
# removing controls as they're already in df2
df1 <- df1[df1$group != "control",]
df <- rbind(df1,df2,df3)

df$dec <- df$dec - df$day
df$time_of_day <- lapply(df$dec, decimal_day_to_hms)
# df$time_of_day <- as_hms(df$time_of_day) # not working
df$time_formatted <- chron(times. = as.character(df$time_of_day),
                           format = "h:m:s")
df$duration_formatted <- chron(times. = as.character(df$feed_duration),
                               format = "h:m:s")
# duration
mean(df$duration_formatted,na.rm = T)
median(df$duration_formatted,na.rm = T)
min(df$duration_formatted,na.rm = T)
max(df$duration_formatted,na.rm = T)

# with/without outliers
test <- df[df$day != 28 & is.na(df$outlier),]
mean(test$duration_formatted,na.rm = T)
median(test$duration_formatted,na.rm = T)
min(test$duration_formatted,na.rm = T)
max(test$duration_formatted,na.rm = T)

# stat specific to day 28
test <- df[df$day == 28,]
mean(test$duration_formatted,na.rm = T)
median(test$duration_formatted,na.rm = T)
min(test$duration_formatted,na.rm = T)
max(test$duration_formatted,na.rm = T)

# time of day
mean(df$time_formatted,na.rm = T)
median(df$time_formatted,na.rm = T)
min(df$time_formatted,na.rm = T)
max(df$time_formatted,na.rm = T)

# with/without outliers
test <- df[df$day != 28 & is.na(df$outlier),]
mean(test$time_formatted,na.rm = T)
median(test$time_formatted,na.rm = T)
min(test$time_formatted,na.rm = T)
max(test$time_formatted,na.rm = T)

# stat specific to day 28
test <- df[df$day == 28,]
mean(test$time_formatted,na.rm = T)
median(test$time_formatted,na.rm = T)
min(test$time_formatted,na.rm = T)
max(test$time_formatted,na.rm = T)
afternoon <- df[df$dec > 0.5,]
sum(afternoon$day == 28) # 13/26 = 50%

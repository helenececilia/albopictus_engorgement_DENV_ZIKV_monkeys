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
library(gammit) # for predict_gamm
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

# Returns the indices of minimum values
which.mins <- function(x, mins=2) {
  head(order(x), mins)
}

# DENV cyno ----

# Old temperature data ----
## Temperature & feeding time data ----
# Temperature ----
temp_denv_cy <- read.csv("../data/DENV_Sylv_Cyno/temperatures.csv")
temp_denv_cy$dec <- temp_denv_cy$Study.Day - as.integer(temp_denv_cy$Study.Day)
temp_denv_cy$time <- case_when(temp_denv_cy$dec == 0 ~ as_hms("00:00:00"),
                               temp_denv_cy$dec == -0.75 ~ as_hms("06:00:00"),
                               temp_denv_cy$dec == -0.5 ~ as_hms("12:00:00"),
                               temp_denv_cy$dec == -0.25 ~ as_hms("18:00:00"),
                               temp_denv_cy$dec == 0.25 ~ as_hms("06:00:00"),
                               temp_denv_cy$dec == 0.5 ~ as_hms("12:00:00"),
                               temp_denv_cy$dec == 0.75 ~ as_hms("18:00:00"))

temp_denv_cy_select <- temp_denv_cy[,c("Study.Day","ID","temp","time")]
temp_denv_cy_select$Day <- as.integer(temp_denv_cy_select$Study.Day)

# Feeding time data ----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_denv_cy <- feeding_time$`DENV Cyno`
feeding_denv_cy <- feeding_denv_cy[,seq(1,5)]
test <- feeding_denv_cy %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                   date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

feed_denv_cy_select <- test[,c("NHP","Day","mid_time")]
colnames(feed_denv_cy_select) <- c("ID","Day","feed_mean")

# Spaces in ID ; matching will not work
feed_denv_cy_select$ID <- gsub(" ","",feed_denv_cy_select$ID, fixed = TRUE)

# Merge and estimate temp at time of feeding 
denv_cy_feed_time_and_temp <- merge(temp_denv_cy_select,feed_denv_cy_select, by = c("ID","Day"))

# We now keep day 0 for a separate model, but beware of negative days which can mess up computations!
# denv_cy_feed_time_and_temp <- denv_cy_feed_time_and_temp[denv_cy_feed_time_and_temp$Day != 0,]
denv_cy_feed_time_and_temp <- denv_cy_feed_time_and_temp[denv_cy_feed_time_and_temp$Study.Day >= 0,]

# denv_cy_feed_time_and_temp$temp_estim_feed_old <- NA

# Fill NAs with mean of known times : for day 0, waiting for data ----
non_na_times <- denv_cy_feed_time_and_temp$feed_mean[!is.na(denv_cy_feed_time_and_temp$feed_mean)]

# To show the existing times were in a narrow window
as_hms(min(non_na_times)) # 8:41
as_hms(max(non_na_times)) # 11:29

non_na_times <- as.character(non_na_times)
mean_time <- mean(times(non_na_times))
mean_time <- as.character(mean_time)
mean_time <- as_hms(mean_time)
denv_cy_feed_time_and_temp$feed_mean[is.na(denv_cy_feed_time_and_temp$feed_mean)] <- mean_time

# Format and save ----
# New method
# More adapted if we have temperatures every 15 mins
# Idea : create columns with the times to do interpolation through group_by / summarize
# but the ordering might not always be the same (depending on which bound the feeding time is closest to)
# so the slope computation might become messy
# ~ 0.21scd  
start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ID", "Day", "temp_estim_feed"))
for(i in unique(denv_cy_feed_time_and_temp$ID)){
  for(d in unique(denv_cy_feed_time_and_temp$Day)){
    test <- denv_cy_feed_time_and_temp[denv_cy_feed_time_and_temp$ID == i &
                                         denv_cy_feed_time_and_temp$Day == d,]
    feed_mean <- unique(test$feed_mean)
    if(!is.na(feed_mean)){
      # find in which window to do linear interpolation
      select <- test[which.mins(abs(test$time - test$feed_mean)),]
      select <- select[order(select$time),] # increasing times
      temp_diff <- select$temp[2] - select$temp[1]
      time_diff <- as.numeric(select$time[2] - select$time[1])
      slope <- temp_diff/time_diff
      temp_estim <- select$temp[1] + 
        slope * as.numeric(feed_mean - select$time[1])
      res <- rbind(res,data.frame(ID = i,
                                  Day = d,
                                  temp_estim_feed = temp_estim))
    }
  }
}
new <- merge(denv_cy_feed_time_and_temp, res, by = c("ID","Day"))
end <- Sys.time()
end-start

# Old method
# Not super optimized yet
# Idea : select parts of the dataframe based on where the feeding time is
# and compute the temperature in a vectorized way
# also avoids going through the NAs for nothing
# ~ 0.56scd
# start <- Sys.time()
# res <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("ID", "Day", "temp_estim_feed_old"))
# for(i in unique(denv_cy_feed_time_and_temp$ID)){
#   for(d in unique(denv_cy_feed_time_and_temp$Day)){
#     df <- denv_cy_feed_time_and_temp[denv_cy_feed_time_and_temp$ID == i &
#                                        denv_cy_feed_time_and_temp$Day == d,]
#     feed_mean <- unique(df$feed_mean)
#     if(!is.na(feed_mean)){
#       if(as_hms("06:00:00") < feed_mean & feed_mean < as_hms("12:00:00")){
#         temp_diff <- df[df$time == as_hms("12:00:00"),]$temp - df[df$time == as_hms("06:00:00"),]$temp
#         time_diff <- as.numeric(as_hms("12:00:00")) - as.numeric(as_hms("06:00:00"))
#         slope <- temp_diff/time_diff
#         temp_estim <- df[df$time == as_hms("06:00:00"),]$temp + 
#           slope * (as.numeric(feed_mean) - as.numeric(as_hms("06:00:00")))
#       }else if(as_hms("12:00:00") < feed_mean & feed_mean < as_hms("18:00:00")){
#         temp_diff <- df[df$time == as_hms("18:00:00"),]$temp - df[df$time == as_hms("12:00:00"),]$temp
#         time_diff <- as.numeric(as_hms("18:00:00")) - as.numeric(as_hms("12:00:00"))
#         slope <- temp_diff/time_diff
#         temp_estim <- df[df$time == as_hms("12:00:00"),]$temp + 
#           slope * (as.numeric(feed_mean) - as.numeric(as_hms("12:00:00")))
#       }else if(as_hms("00:00:00") < feed_mean & feed_mean < as_hms("06:00:00")){
#         temp_diff <- df[df$time == as_hms("06:00:00"),]$temp - df[df$time == as_hms("00:00:00"),]$temp
#         time_diff <- as.numeric(as_hms("06:00:00")) - as.numeric(as_hms("00:00:00"))
#         slope <- temp_diff/time_diff
#         temp_estim <- df[df$time == as_hms("00:00:00"),]$temp + 
#           slope * (as.numeric(feed_mean) - as.numeric(as_hms("00:00:00")))
#       }else if(as_hms("18:00:00") < feed_mean){
#         # need to do the computation with 00:00:00 of the next day
#         print("Feeding after 18:00")
#         print("ID")
#         print(i)
#         print("Day")
#         print(d)
#       }
#       res <- rbind(res,data.frame(ID = i,
#                                   Day = d,
#                                   temp_estim_feed_old = temp_estim))
#     }
#   }
# }
# old <- merge(denv_cy_feed_time_and_temp, res, by = c("ID","Day"))
# end <- Sys.time()
# end-start

write.csv(new, "../output/mosq_feeding_behaviour/DENV_Cyno_temperature_estimate_at_feeding.csv", row.names = F)

## Converting feeding times into day and its decimal ----
denv_cy_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/DENV_Cyno_temperature_estimate_at_feeding.csv")

denv_cy_feed_time_and_temp$chron_time <- chron(times. = as.character(denv_cy_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

denv_cy_feed_time_and_temp$dec <- denv_cy_feed_time_and_temp$Day + hms_to_decimal_day(denv_cy_feed_time_and_temp$chron_time)

## Some plots ----
ggplot(denv_cy_feed_time_and_temp) + 
  geom_point(aes(x = Study.Day, y = temp),
             color = "blue", alpha = 0.5) +
  geom_line(aes(x = Study.Day, y = temp),
            color = "blue", alpha = 0.5) +
  geom_point(aes(x = dec, y = temp_estim_feed),
             color = "orange") + 
  facet_wrap(~ ID, scales = "free_y")

## Including feeding behaviour ----
colnames(denv_cy_feed_time_and_temp)[2] <- "day"

select_temp_cyno <- unique(denv_cy_feed_time_and_temp[,c("ID","day","temp_estim_feed","dec")])

df <- read.csv("../data/DENV_Sylv_Cyno/Trade-Off_Project_Cyno_Inf_Sylvatic_DENV-2_2021_Comprehensive_data_27Jan23_HC.csv", sep = "\t", dec = ".") # , dec = "," for previous file version (17Nov BM)

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
day_0$total <- ifelse(day_0$nb_fed_day0 > 1, 10, 1)
day_0_with_temp <- merge(day_0, select_temp_cyno, by = c("ID","day"))
write.csv(day_0_with_temp, "../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_day0.csv",
          row.names = F)

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
others <- others[complete.cases(others),] 
denv_cy <- merge(others, select_temp_cyno, by = c("ID","day"))

# 10 is the default, occasionnally more
denv_cy$total <- ifelse(denv_cy$nb_fed <= 10, 10, denv_cy$nb_fed)

write.csv(denv_cy, "../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour.csv", row.names = F)


# ZIKV squirrel ------
## Temperature & feeding time data ----
# Temperature ----
temp_zikv <- read.csv("../data/ZIKV_Sylv_Squirrel/temperatures.csv")

temp_zikv$dec <- temp_zikv$Study.Day - as.integer(temp_zikv$Study.Day)
temp_zikv$time <- case_when(temp_zikv$dec == 0 ~ as_hms("00:00:00"),
                            temp_zikv$dec == -0.75 ~ as_hms("06:00:00"),
                            temp_zikv$dec == -0.5 ~ as_hms("12:00:00"),
                            temp_zikv$dec == -0.25 ~ as_hms("18:00:00"),
                            temp_zikv$dec == 0.25 ~ as_hms("06:00:00"),
                            temp_zikv$dec == 0.5 ~ as_hms("12:00:00"),
                            temp_zikv$dec == 0.75 ~ as_hms("18:00:00"))

temp_zikv_select <- temp_zikv[,c("Study.Day","ID","temp","time")]
temp_zikv_select$Day <- as.integer(temp_zikv_select$Study.Day)

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

feed_zikv_select <- test[,c("NHP","Day","mid_time","mosq_type")]
colnames(feed_zikv_select) <- c("ID","Day","feed_mean","mosq_type")

# Merge and estimate temp at time of feeding
zikv_feed_time_and_temp <- merge(temp_zikv_select,feed_zikv_select, by = c("ID","Day"))
# We keep day 0 now, but remove < 0
zikv_feed_time_and_temp <- zikv_feed_time_and_temp[zikv_feed_time_and_temp$Study.Day >= 0,]
zikv_feed_time_and_temp <- zikv_feed_time_and_temp[complete.cases(zikv_feed_time_and_temp),]

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("ID", "Day", "temp_estim_feed","mosq_type"))
for(i in unique(zikv_feed_time_and_temp$ID)){
  for(d in unique(zikv_feed_time_and_temp$Day)){
    test <- zikv_feed_time_and_temp[zikv_feed_time_and_temp$ID == i &
                                      zikv_feed_time_and_temp$Day == d,]
    for(t in unique(test$mosq_type)){ # for day 0, 2 different types, each with their own computation
      df <- test[test$mosq_type == t,]
      feed_mean <- unique(df$feed_mean)
      # if(length(feed_mean)>1){browser()}
      if(!is.na(feed_mean)){
        # find in which window to do linear interpolation
        select <- df[which.mins(abs(df$time - df$feed_mean)),]
        select <- select[order(select$time),] # increasing times
        temp_diff <- select$temp[2] - select$temp[1]
        time_diff <- as.numeric(select$time[2] - select$time[1])
        slope <- temp_diff/time_diff
        temp_estim <- select$temp[1] + 
          slope * as.numeric(feed_mean - select$time[1])
        res <- rbind(res,data.frame(ID = i,
                                    Day = d,
                                    temp_estim_feed = temp_estim,
                                    mosq_type = t))
      }
    }
  }
}
new <- merge(zikv_feed_time_and_temp, res, by = c("ID","Day","mosq_type"))
end <- Sys.time()
end-start

write.csv(new, "../output/mosq_feeding_behaviour/ZIKV_Squirrel_temperature_estimate_at_feeding.csv",
          row.names = F)

## Converting feeding times into day and its decimal ----
zikv_sq_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/ZIKV_Squirrel_temperature_estimate_at_feeding.csv")

zikv_sq_feed_time_and_temp$chron_time <- chron(times. = as.character(zikv_sq_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

zikv_sq_feed_time_and_temp$dec <- zikv_sq_feed_time_and_temp$Day + hms_to_decimal_day(zikv_sq_feed_time_and_temp$chron_time)

## Some plots ----
ggplot(zikv_sq_feed_time_and_temp) + 
  geom_point(aes(x = Study.Day, y = temp),
             color = "blue", alpha = 0.5) +
  geom_line(aes(x = Study.Day, y = temp),
            color = "blue", alpha = 0.5) +
  geom_point(aes(x = dec, y = temp_estim_feed),
             color = "orange") + 
  facet_wrap(~ ID, scales = "free_y")


# DENV squirrel -----
## Temperature & feeding time data ----
# Temperature ----
temp_denv_sq <- read.csv("../data/DENV_Sylv_Squirrel/temperatures.csv")
temp_denv_sq$dec <- temp_denv_sq$Study.Day - as.integer(temp_denv_sq$Study.Day)
temp_denv_sq$time <- case_when(temp_denv_sq$dec == 0 ~ as_hms("00:00:00"),
                               temp_denv_sq$dec == -0.75 ~ as_hms("06:00:00"),
                               temp_denv_sq$dec == -0.5 ~ as_hms("12:00:00"),
                               temp_denv_sq$dec == -0.25 ~ as_hms("18:00:00"),
                               temp_denv_sq$dec == 0.25 ~ as_hms("06:00:00"),
                               temp_denv_sq$dec == 0.5 ~ as_hms("12:00:00"),
                               temp_denv_sq$dec == 0.75 ~ as_hms("18:00:00"))

temp_denv_sq_select <- temp_denv_sq[,c("Study.Day","ID","temp","time")]
temp_denv_sq_select$Day <- as.integer(temp_denv_sq_select$Study.Day)

# Feeding time data ----
feeding_time <- multiplesheets("../data/Trade-Off Project. Squirrel Monkeys AND Cynomologous Macaques. Mosquito feeding experiment times_01_Feb_2023_HC.xlsx")
feeding_denv_sq <- feeding_time$`DENV Saimiri`
feeding_denv_sq$Notes <- NULL
test <- feeding_denv_sq %>% mutate(date1 = ymd_hms(`Mosquito feed start`),
                                   date2 = ymd_hms(`Mosquito feed stop`))
test$hour1<- as_hms(test$date1)
test$hour2 = as_hms(test$date2)
test$diff = test$hour2 - test$hour1 
test$mid_time <- as_hms(test$hour1 + test$diff/2)

feed_denv_sq_select <- test[,c("NHP","Day","mid_time","mosq_type")]
colnames(feed_denv_sq_select) <- c("ID","Day","feed_mean","mosq_type")

# Merge and estimate temp at time of feeding
denv_sq_feed_time_and_temp <- merge(temp_denv_sq_select,feed_denv_sq_select, by = c("ID","Day"))
# We keep day 0 now, but remove < 0
denv_sq_feed_time_and_temp <- denv_sq_feed_time_and_temp[denv_sq_feed_time_and_temp$Study.Day >= 0,]
denv_sq_feed_time_and_temp <- denv_sq_feed_time_and_temp[complete.cases(denv_sq_feed_time_and_temp),]

start <- Sys.time()
res <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("ID", "Day", "temp_estim_feed","mosq_type"))
for(i in unique(denv_sq_feed_time_and_temp$ID)){
  for(d in unique(denv_sq_feed_time_and_temp$Day)){
    test <- denv_sq_feed_time_and_temp[denv_sq_feed_time_and_temp$ID == i &
                                         denv_sq_feed_time_and_temp$Day == d,]
    for(t in unique(test$mosq_type)){ # for day 0, 2 different types, each with their own computation
      df <- test[test$mosq_type == t,]
      feed_mean <- unique(df$feed_mean)
      if(!is.na(feed_mean)){
        # find in which window to do linear interpolation
        select <- df[which.mins(abs(df$time - df$feed_mean)),]
        select <- select[order(select$time),] # increasing times
        temp_diff <- select$temp[2] - select$temp[1]
        time_diff <- as.numeric(select$time[2] - select$time[1])
        slope <- temp_diff/time_diff
        temp_estim <- select$temp[1] + 
          slope * as.numeric(feed_mean - select$time[1])
        res <- rbind(res,data.frame(ID = i,
                                    Day = d,
                                    temp_estim_feed = temp_estim,
                                    mosq_type = t))
      }
    }
  }
}
new <- merge(denv_sq_feed_time_and_temp, res, by = c("ID","Day","mosq_type"))
end <- Sys.time()
end-start

write.csv(new, "../output/mosq_feeding_behaviour/DENV_Squirrel_temperature_estimate_at_feeding.csv",
          row.names = F)

## Converting feeding times into day and its decimal ----
denv_sq_feed_time_and_temp <- read.csv("../output/mosq_feeding_behaviour/DENV_Squirrel_temperature_estimate_at_feeding.csv")

denv_sq_feed_time_and_temp$chron_time <- chron(times. = as.character(denv_sq_feed_time_and_temp$feed_mean),
                                               format = "h:m:s")

denv_sq_feed_time_and_temp$dec <- denv_sq_feed_time_and_temp$Day + hms_to_decimal_day(denv_sq_feed_time_and_temp$chron_time)

## Some plots ----
ggplot(denv_sq_feed_time_and_temp) + 
  geom_point(aes(x = Study.Day, y = temp),
             color = "blue", alpha = 0.5) +
  geom_line(aes(x = Study.Day, y = temp),
            color = "blue", alpha = 0.5) +
  geom_point(aes(x = dec, y = temp_estim_feed),
             color = "orange") + 
  facet_wrap(~ ID, scales = "free_y")

# Including feeding behaviour ----
# Aggregated controls squirrels in zikv and denv datasets
colnames(zikv_sq_feed_time_and_temp)[2] <- "day"
colnames(denv_sq_feed_time_and_temp)[2] <- "day"

temp_squirrel <- rbind(zikv_sq_feed_time_and_temp, denv_sq_feed_time_and_temp)
select_temp_squirrel <- unique(temp_squirrel[,c("ID","day","mosq_type","temp_estim_feed","dec")])

df <- read.csv("../data/ZIKV_Sylv_Squirrel/Trade-Off_Project_Saimiri_Inf_Sylvatic_ZIKV_2021_Comprehensive_data_25Jan2023_BM.csv", dec = ".", sep = "\t")
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
# write.csv(day_0_with_temp, "../output/mosq_feeding_behaviour/DENV_Cyno_data_feeding_behaviour_day0.csv",
#           row.names = F)

# Dataset other days
others <- df[df$day > 0,]
others$nb_fed_day0 <- NULL
others <- others[complete.cases(others),] 
zikv_inf_others <- merge(others, select_temp_squirrel, by = c("ID","day"))


df <- read.csv("../data/DENV_Sylv_Squirrel/Trade-Off_Project_Saimiri_Inf_Sylvatic_DENV-2_2021_Comprehensive_data_19Jan23_BM.csv", dec = ".", sep = "\t")
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

controls_zikv_sq <- read.csv("../data/ZIKV_Sylv_Squirrel/Trade-Off_Project_Saimiri_Inf_Sylvatic_ZIKV_2021_Comprehensive_data_25Jan2023_BM.csv", dec = ".", sep = "\t")
controls_zikv_sq <- controls_zikv_sq[controls_zikv_sq$Final.Treatment == "Control",]

controls_denv_sq <- read.csv("../data/DENV_Sylv_Squirrel/Trade-Off_Project_Saimiri_Inf_Sylvatic_DENV-2_2021_Comprehensive_data_19Jan23_BM.csv", dec = ".", sep = "\t")
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

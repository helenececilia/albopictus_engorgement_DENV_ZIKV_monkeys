## ---------------------------
## Author: Helene Cecilia
##
## Date Created: 2023-04-05

rm(list=ls())

## Loading Packages  ------------------
library(ggdist)
library(tidyverse)
library(patchwork)
library(Matching) # for ks.boot
library(scales)
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
library(janitor) # for clean_names
library(multcomp) # for glht
library(gghalves) # for half_boxplot
library(cowplot)
library(magick)
library(ggtext) # for element_markdown


## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Utility functions ---

hms_to_decimal_day <-function(time_hms){
  h <- hours(time_hms)
  m <- minutes(time_hms)
  s <- seconds(time_hms)
  sum_minutes <- h*60 + m + s/60
  dec <- sum_minutes / (24*60)
  return(dec)
}

# Engorgement data after day 0 ----
zikv_sq_data <- read.csv("../data/ZIKV_Squirrel_data_engorgement.csv")
zikv_sq_data$species <- "Squirrel"
zikv_sq_data$virus <- "Zika"
zikv_sq_data$virus[zikv_sq_data$group == "control"] <- "none"

zikv_cy_data <- read.csv("../data/ZIKV_Cyno_data_engorgement.csv")
zikv_cy_data$species <- "Cyno"
zikv_cy_data$virus <- "Zika"

denv_sq_data <- read.csv("../data/DENV_Squirrel_data_engorgement.csv")
denv_sq_data$species <- "Squirrel"
denv_sq_data$virus <- "Dengue"
denv_sq_data$virus[denv_sq_data$group == "control"] <- "none"

denv_cy_data <- read.csv("../data/DENV_Cyno_data_engorgement.csv")
denv_cy_data$species <- "Cyno"
denv_cy_data$virus <- "Dengue"
denv_cy_data$virus[denv_cy_data$group == "control"] <- "none"

my_data <- rbind(zikv_sq_data, zikv_cy_data, denv_sq_data, denv_cy_data)
my_data <- unique(my_data)

my_data$time_feed <- my_data$dec - my_data$day

my_data$sex <- as.factor(my_data$sex)
my_data <- within(my_data, sex <- relevel(sex, ref = "M"))
my_data$species <- as.factor(my_data$species)
my_data <- within(my_data, species <- relevel(species, ref = "Squirrel"))
my_data$virus <- as.factor(my_data$virus)
my_data <- within(my_data, virus <- relevel(virus, ref = "Dengue"))

my_data$int_duration <- minutes(chron(times. = my_data$feed_duration,
                                      format = "h:m:s"))
eng1 <- my_data # 331 obs

# Temperature data ----
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

df1$Study.Day <- df1$Day + hms_to_decimal_day(df1$chron_time)
df2$Study.Day <- df2$Day + hms_to_decimal_day(df2$chron_time)
df3$Study.Day <- df3$Day + hms_to_decimal_day(df3$chron_time)
df4$Study.Day <- df4$Day + hms_to_decimal_day(df4$chron_time)

# max(c(df1$temp,df2$temp,df3$temp,df4$temp)) # 40.65

df1$group <- "control"
df1$group[df1$Final.Treatment == "1 Mosquito"] <- "infected"
df1$group[df1$Final.Treatment == "10 Mosquitos"] <- "infected"
df2$group <- "control"
df2$group[df2$Final.Treatment == "15 Mosquitos"] <- "infected"
df3$group <- "control"
df3$group[df3$Final.Treatment == "15 Mosquitos"] <- "infected"
df4$group <- "control"
df4$group[df4$Final.Treatment == "15 Mosquitos"] <- "infected"

# Data selection and formatting ----
# Only show temperature of days when feeding happened
colnames(eng1)[2] <- "Day"
df1 <- df1[df1$Day %in% unique(eng1[eng1$species == "Cyno",]$Day),]
df2 <- df2[df2$Day %in% unique(eng1[eng1$species == "Squirrel" & 
                                      eng1$virus %in% c("Dengue","none"),]$Day),]
df3 <- df3[df3$Day %in% unique(eng1[eng1$species == "Squirrel" & 
                                      eng1$virus %in% c("Zika","none"),]$Day),]
df4 <- df4[df4$Day %in% unique(eng1[eng1$species == "Cyno" & 
                                      eng1$virus %in% c("Zika","none"),]$Day),]

df1_cut <- merge(df1,eng1[eng1$species == "Cyno",
                          c("ID","Day","dec")], by = c("ID","Day"))

df2_cut <- merge(df2,eng1[eng1$species == "Squirrel" & 
                            eng1$virus %in% c("Dengue","none"),
                          c("ID","Day","dec")], by = c("ID","Day"))

df3_cut <- merge(df3,eng1[eng1$species == "Squirrel" & 
                            eng1$virus %in% c("Zika","none"),
                          c("ID","Day","dec")], by = c("ID","Day"))

df4_cut <- merge(df4,eng1[eng1$species == "Cyno" & 
                            eng1$virus %in% c("Zika","none"),
                          c("ID","Day","dec")], by = c("ID","Day"))

# Remove temperature dynamics of day 28 after feeding to not see T° going down because of euthanasia
df1_cut$temp[df1_cut$Day == 28 & df1_cut$Study.Day > df1_cut$dec + 0.01] <- NA
df2_cut$temp[df2_cut$Day == 28 & df2_cut$Study.Day > df2_cut$dec + 0.01] <- NA
df3_cut$temp[df3_cut$Day == 28 & df3_cut$Study.Day > df3_cut$dec + 0.01] <- NA
df4_cut$temp[df4_cut$Day == 28 & df4_cut$Study.Day > df4_cut$dec + 0.01] <- NA

# Remove day 0 
df1_cut <- df1_cut[df1_cut$Day != 0,]
df2_cut <- df2_cut[df2_cut$Day != 0,]
df3_cut <- df3_cut[df3_cut$Day != 0,]
df4_cut <- df4_cut[df4_cut$Day != 0,]

df1_cut$monkey_status <- df1_cut$group
df1_cut$monkey_status[df1_cut$group == "control"] <- "Control"
df1_cut$monkey_status[df1_cut$group != "control"] <- "Dengue-2 virus"

df2_cut$monkey_status <- df2_cut$group
df2_cut$monkey_status[df2_cut$group == "control"] <- "Control"
df2_cut$monkey_status[df2_cut$group != "control"] <- "Dengue-2 virus"

df3_cut$monkey_status <- df3_cut$group
df3_cut$monkey_status[df3_cut$group == "control"] <- "Control"
df3_cut$monkey_status[df3_cut$group != "control"] <- "Zika virus"

df4_cut$monkey_status <- df4_cut$group
df4_cut$monkey_status[df4_cut$group == "control"] <- "Control"
df4_cut$monkey_status[df4_cut$group != "control"] <- "Zika virus"

eng1$monkey_status <- eng1$virus
eng1$monkey_status <- as.character(eng1$monkey_status)
eng1$monkey_status[eng1$monkey_status == "none"] <- "Control"
eng1$monkey_status[eng1$monkey_status == "Dengue"] <- "Dengue-2 virus"
eng1$monkey_status[eng1$monkey_status == "Zika"] <- "Zika virus"

# to plot control points on top (not working?)
eng1$monkey_status <- factor(eng1$monkey_status,
                             levels = c("Dengue-2 virus", "Zika virus", "Control"))
df1_cut$monkey_status <- factor(df1_cut$monkey_status,
                                levels = c("Dengue-2 virus", "Control"))
df2_cut$monkey_status <- factor(df2_cut$monkey_status,
                                levels = c("Dengue-2 virus", "Control"))
df3_cut$monkey_status <- factor(df3_cut$monkey_status,
                                levels = c("Zika virus", "Control"))
df4_cut$monkey_status <- factor(df4_cut$monkey_status,
                                levels = c("Zika virus", "Control"))

# Temperature dynamics : days superposed -----
pd1 <- ggplot() + geom_line(data = df1_cut,
                            aes(x = Study.Day - Day, y = temp,
                                color = monkey_status, group = interaction(ID,Day)),
                            alpha = 0.2) +
  geom_point(data = eng1[eng1$species == "Cyno" & 
                           eng1$virus %in% c("Dengue","none") &
                           eng1$Day == 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 shape = monkey_status),
             fill = "black",
             color = "red",
             size = 3, alpha = 0.75) +
  geom_point(data = eng1[eng1$species == "Cyno" & 
                           eng1$virus %in% c("Dengue","none") &
                           eng1$Day != 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 fill = monkey_status,
                 shape = monkey_status),
             color = "#414341",
             size = 4, alpha = 0.4) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#C2c3c9")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  labs(x = "", y = expression("Body temperature " ( degree*C)),
       color = "Monkey infection status", fill = "Monkey infection status", shape = "Monkey infection status") +
  coord_cartesian(xlim = c(0,1), ylim = c(34,40.7)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.75)),
         fill = guide_legend(override.aes = list(size = 3.5))) +
  ggtitle("A") +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = c("00:00","6AM","12:00","6PM","00:00"),
                     expand = expansion(add = c(0.01,0.01))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 29),
        axis.title.y = element_text(size = 29,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 27),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28), 
        plot.title = element_text(size = 27),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(0,26,0,0))

pd2 <- ggplot() + geom_line(data = df2_cut,
                            aes(x = Study.Day - Day, y = temp,
                                color = monkey_status, group = interaction(ID,Day)),
                            alpha = 0.2) +
  geom_point(data = eng1[eng1$species == "Squirrel" & 
                           eng1$virus %in% c("Dengue","none") & 
                           eng1$Day == 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 shape = monkey_status),
             fill = "black",
             color = "red",
             size = 3, alpha = 0.75) +
  geom_point(data = eng1[eng1$species == "Squirrel" & 
                           eng1$virus %in% c("Dengue","none") & 
                           eng1$Day != 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 fill = monkey_status,
                 shape = monkey_status),
             color = "#414341",
             size = 4, alpha = 0.4) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#C2c3c9")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  labs(x = "", y = "",
       color = "", fill = "", shape = "") +
  coord_cartesian(xlim = c(0,1), ylim = c(34,40.7)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.75)),
         fill = guide_legend(override.aes = list(size = 3.5))) +
  ggtitle("B") +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = c("00:00","6AM","12:00","6PM","00:00"),
                     expand = expansion(add = c(0.01,0.01))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 29),
        axis.title.y = element_text(size = 29,
                                    margin = margin(r = 1)),
        axis.text = element_text(size = 27),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 27),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(0,22,0,0))

pd3 <- ggplot() + geom_line(data = df3_cut,
                            aes(x = Study.Day - Day, y = temp,
                                color = monkey_status, group = interaction(ID,Day)),
                            alpha = 0.2) +
  geom_point(data = eng1[eng1$species == "Squirrel" & 
                           eng1$virus %in% c("Zika","none") & 
                           eng1$Day == 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 shape = monkey_status),
             fill = "black",
             color = "red",
             size = 3, alpha = 0.75) +
  geom_point(data = eng1[eng1$species == "Squirrel" & 
                           eng1$virus %in% c("Zika","none") & 
                           eng1$Day != 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 fill = monkey_status,
                 shape = monkey_status),
             color = "#414341",
             size = 4, alpha = 0.4) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#C2c3c9")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  labs(x = "Time of day", y = "",
       color = "", fill = "", shape = "") +
  coord_cartesian(xlim = c(0,1), ylim = c(34,40.7)) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = c("00:00","6AM","12:00","6PM","00:00"),
                     expand = expansion(add = c(0.01,0.01))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.75)),
         fill = guide_legend(override.aes = list(size = 3.5))) +
  ggtitle("D") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 29,
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 29,
                                    margin = margin(r = 1)),
        axis.text = element_text(size = 27),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 27),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(0,22,0,0))

pd4 <- ggplot() + geom_line(data = df4_cut,
                            aes(x = Study.Day - Day, y = temp,
                                color = monkey_status, group = interaction(ID,Day)),
                            alpha = 0.2) +
  geom_point(data = eng1[eng1$species == "Cyno" & 
                           eng1$virus %in% c("Zika","none") & 
                           eng1$Day == 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 shape = monkey_status),
             fill = "black",
             color = "red",
             size = 3, alpha = 0.75) +
  geom_point(data = eng1[eng1$species == "Cyno" & 
                           eng1$virus %in% c("Zika","none") & 
                           eng1$Day != 28,],
             aes(x = time_feed, y = temp_estim_feed,
                 fill = monkey_status,
                 shape = monkey_status),
             color = "#414341",
             size = 4, alpha = 0.4) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#C2c3c9")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  labs(x = "Time of day", y = expression("Body temperature " ( degree*C)),
       color = "", fill = "", shape = "") +
  coord_cartesian(xlim = c(0,1), ylim = c(34,40.7)) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = c("00:00","6AM","12:00","6PM","00:00"),
                     expand = expansion(add = c(0.01,0.01))) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.75)),
         fill = guide_legend(override.aes = list(size = 3.5))) +
  ggtitle("C") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 29,
                                    margin = margin(t = 15)),
        axis.title.y = element_text(size = 29,
                                    margin = margin(r = 15)),
        axis.text = element_text(size = 27),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 27),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(0,26,0,0))

img_cyno <- magick::image_read("../output/figures/outlines/macaque_outline.png")
img_squirrel <- magick::image_read("../output/figures/outlines/squirrel_outline.png")

p1 <- ggdraw() + 
  draw_plot(pd1) +
  draw_image(image = img_cyno, 
             x = 0.15, y = 0.71, scale = 0.15,
             valign = 0, halign = 0)

p2 <- ggdraw() + 
  draw_plot(pd2) +
  draw_image(image = img_squirrel, 
             x = 0.12, y = 0.71, scale = 0.15,
             valign = 0, halign = 0)

p3 <- ggdraw() + 
  draw_plot(pd3) +
  draw_image(image = img_squirrel, 
             x = 0.12, y = 0.78, scale = 0.15,
             valign = 0, halign = 0)

p4 <- ggdraw() + 
  draw_plot(pd4) +
  draw_image(image = img_cyno, 
             x = 0.15, y = 0.78, scale = 0.15,
             valign = 0, halign = 0)


pd <- (p1|p2)/(p4|p3)

png(filename = "../output/figures/Figure_S1.png",
    width = 1950, height = 1300)
print(pd)
dev.off()

# Temperature and engorgement - without day 28 -----
eng1 <- eng1[eng1$Day != 28,] 

eng1$prop <- eng1$nb_fed/eng1$total

ptt1 <- ggplot(eng1[eng1$species == "Cyno" & eng1$virus %in% c("none","Dengue"),]) +
  geom_point(aes(y = prop,
                 x = temp_estim_feed,
                 color = monkey_status,
                 fill = monkey_status,
                 shape = monkey_status),
             size = 5,
             alpha = 0.65,
             stroke = 1.2) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "grey")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#414341")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(min(eng1$temp_estim_feed), max(eng1$temp_estim_feed)),
                  clip = "off") + # to allow text outside plot area
  scale_x_continuous(breaks = c(36,37,38,39,40)) +
  labs(x = "",
       y = "Proportion of engorged mosquitoes",
       fill = bquote(bold("Monkey infection status")),
       color = bquote(bold("Monkey infection status")),
       shape = bquote(bold("Monkey infection status"))) +
  ggtitle("A") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(alpha = 1, size = 5)),
         fill = guide_legend(order = 1),
         color = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 1)),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        plot.margin = margin(2,10,0,10),
        plot.title = element_text(size = 27))

ptt2 <- ggplot(eng1[eng1$species == "Squirrel" & eng1$virus %in% c("none","Dengue"),]) +
  geom_point(aes(y = prop,
                 x = temp_estim_feed,
                 color = monkey_status,
                 fill = monkey_status,
                 shape = monkey_status),
             size = 5,
             alpha = 0.65,
             stroke = 1.2) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "grey")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#414341")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(min(eng1$temp_estim_feed), max(eng1$temp_estim_feed))) +
  scale_x_continuous(breaks = c(36,37,38,39,40)) +
  labs(x = "",
       y = "",
       fill = "", color = "", shape = "") +
  ggtitle("B") +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 1)),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 1)),
        axis.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        plot.margin = margin(2,10,0,10),
        plot.title = element_text(size = 27))

ptt3 <- ggplot(eng1[eng1$species == "Squirrel" & eng1$virus %in% c("none","Zika"),]) +
  geom_point(aes(y = prop,
                 x = temp_estim_feed,
                 color = monkey_status,
                 fill = monkey_status,
                 shape = monkey_status),
             size = 5,
             alpha = 0.65,
             stroke = 1.2) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "grey")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#414341")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(min(eng1$temp_estim_feed), max(eng1$temp_estim_feed))) +
  scale_x_continuous(breaks = c(36,37,38,39,40)) +
  labs(x = expression("Body temperature " ( degree*C)),
       y = "",
       fill = "", color = "", shape = "") +
  ggtitle("D") +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 10)),
        axis.title.y = element_text(size = 30,
                                    margin = margin(r = 1)),
        axis.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.margin = margin(2,10,0,10),
        plot.title = element_text(size = 27))

ptt4 <- ggplot(eng1[eng1$species == "Cyno" & eng1$virus %in% c("none","Zika"),]) +
  geom_point(aes(y = prop,
                 x = temp_estim_feed,
                 color = monkey_status,
                 fill = monkey_status,
                 shape = monkey_status),
             size = 5,
             alpha = 0.65,
             stroke = 1.2) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "grey")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#414341")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(min(eng1$temp_estim_feed), max(eng1$temp_estim_feed))) +
  scale_x_continuous(breaks = c(36,37,38,39,40)) +
  labs(x = expression("Body temperature " ( degree*C)),
       y = "Proportion of engorged mosquitoes",
       fill = "", color = "", shape = "") +
  ggtitle("C") +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30,
                                    margin = margin(t = 10)),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 15)),
        axis.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.margin = margin(2,10,0,10),
        plot.title = element_text(size = 27))

img_cyno <- magick::image_read("../output/figures/outlines/macaque_outline.png")
img_squirrel <- magick::image_read("../output/figures/outlines/squirrel_outline.png")

p1 <- ggdraw() + 
  draw_plot(ptt1) +
  draw_image(image = img_cyno, 
             x = 0.85, y = 0.12, scale = 0.15,
             valign = 0, halign = 0)

p2 <- ggdraw() + 
  draw_plot(ptt2) +
  draw_image(image = img_squirrel, 
             x = 0.85, y = 0.12, scale = 0.15,
             valign = 0, halign = 0)

p3 <- ggdraw() + 
  draw_plot(ptt3) +
  draw_image(image = img_squirrel, 
             x = 0.85, y = 0.22, scale = 0.15,
             valign = 0, halign = 0)

p4 <- ggdraw() + 
  draw_plot(ptt4) +
  draw_image(image = img_cyno, 
             x = 0.85, y = 0.22, scale = 0.15,
             valign = 0, halign = 0)

p <- (p1/p4)|(p2/p3)

png(filename = "../output/figures/Figure_4.png",
    width = 1950, height = 1300)
print(p)
dev.off()


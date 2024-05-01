## ---------------------------
## Author: Helene Cecilia
##
## Date Created: 2023-02-10

rm(list=ls())

## Loading Packages  ------------------
library(tidyverse)
library(BAS)
library(patchwork)
library(cowplot) # for draw_image
library(magick) # for draw_image / note : installed through command line, error in Rstudio
library(ggpubr) # for stat_compare_means, geom_bracket
library(car) # for vif and anova
library(performance)
library(scales) # for trans_breaks, alpha function
library(broom) # for tidy
library(DHARMa)
library(effects)
library(glmmTMB)
library(data.table) # for %like%
library(MuMIn) # for AICc to work on glmmTMB objects
library(multcomp) # for glht
library(gghalves)
library(extraDistr) # for beta-binomial (pbbinom, qbbinom etc.)
library(emdbook) # for beta-binomial (rbetabinom, dbetabinom)
library(piecewiseSEM) # for dagify
library(ggdag) # to plot dag
library(chron) # for minutes
library(lsr) # for etaSquared
library(buildmer) # for buildglmmTMB / remotes::install_github('cvoeten/buildmer')
library(ggtext) # for element_markdown
library(ggeffects) # for geom_errorbar

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# The Data -----
# Include cyno high exposure only
denv_cy <- read.csv("../data/Day0_DENV_Cyno_data_engorgement.csv")
denv_cy <- denv_cy[denv_cy$total != 1,]
denv_cy$mosq_type <- ifelse(denv_cy$group == "infected","infecting","control")
denv_cy$species <- "Cyno"
denv_cy$virus_true <- ifelse(denv_cy$mosq_type == "infecting","Dengue virus","")
denv_cy$virus_label <- denv_cy$virus_true

denv_sq <- read.csv("../data/Day0_DENV_Squirrel_data_engorgement.csv")
denv_sq$species <- "Squirrel"
# two categories : "virus_true" used to quantify virus effect 
# (concurrent feeding mosquitoes should not 
# be aggregated with infected mosquitoes for a same virus)
# "virus_label" used for labeling the groups : 
# concurrent feeding associated with a virus to differenciate
denv_sq$virus_true <- ifelse(denv_sq$mosq_type == "infecting","Dengue virus","")
denv_sq$virus_label <- ifelse(denv_sq$mosq_type == "control","","Dengue virus")

zikv_sq <- read.csv("../data/Day0_ZIKV_Squirrel_data_engorgement.csv")
zikv_sq$species <- "Squirrel"
zikv_sq$virus_true <- ifelse(zikv_sq$mosq_type == "infecting","Zika virus","")
zikv_sq$virus_label <- ifelse(zikv_sq$mosq_type == "control","","Zika virus")

zikv_cy <- read.csv("../data/Day0_ZIKV_Cyno_data_engorgement.csv")
zikv_cy$species <- "Cyno"
zikv_cy$mosq_type <- ifelse(zikv_cy$group == "infected","infecting","control")
zikv_cy$virus_true <- ifelse(zikv_cy$mosq_type == "infecting","Zika virus","")
zikv_cy$virus_label <- ifelse(zikv_cy$mosq_type == "control","","Zika virus")

df_virus <- rbind(zikv_cy,zikv_sq,denv_sq,denv_cy)
df_virus$mosq_type_agg <- interaction(df_virus$mosq_type,df_virus$virus_label)

df_virus$species <- as.factor(df_virus$species)
df_virus <- within(df_virus, species <- relevel(species, ref = "Squirrel"))

df_virus$mosq_type_agg <- as.factor(df_virus$mosq_type_agg)
df_virus <- within(df_virus, mosq_type_agg <- relevel(mosq_type_agg, ref = "infecting.Dengue virus"))

df_virus$sex <- as.factor(df_virus$sex)
df_virus <- within(df_virus, sex <- relevel(sex, ref = "F"))

df_virus$weight <- as.numeric(df_virus$weight)

df_virus <- unique(df_virus) # remove duplicate controls

df_virus$mosq_nhp <- interaction(df_virus$species,df_virus$mosq_type_agg)
df_virus$mosq_nhp <- gsub(" ","_",df_virus$mosq_nhp) # no space for posthoc test to work
df_virus$mosq_nhp <- as.factor(df_virus$mosq_nhp)
df_virus <- within(df_virus, mosq_nhp <- relevel(mosq_nhp, ref = "Squirrel.infecting.Dengue_virus"))

df_virus$time_feed <- df_virus$dec - df_virus$day

feed_duration <- str_split(df_virus$feed_duration,":")
vec_duration <- sapply(feed_duration,"[[",2)
df_virus$int_duration <- as.integer(vec_duration)

# a <- df_virus$int_duration[df_virus$species == "Cyno"]
# summary(a)
# sd(a)/sqrt(length(a))
# 
# b <- df_virus$int_duration[df_virus$species == "Squirrel" & df_virus$mosq_type != "cofeeding"]
# summary(b)
# sd(b)/sqrt(length(b))
# 
# c <- df_virus$int_duration[df_virus$mosq_type == "cofeeding"]
# summary(c)
# sd(c)/sqrt(length(c))

test <- df_virus
test$group_agg <- interaction(df_virus$species,df_virus$virus_true)
test$group_agg <- gsub(" ","_",test$group_agg) # no space for posthoc test to work

test <- test[,c("ID","nb_fed_day0","total","virus_true","species","group_agg","int_duration")]
# Add the unique concurrent feeding batch on cynos
test <- rbind(test,list("FR423A",45,50,"","Cyno","Cyno.",NA)) # we don't know the duration for this feeding 
# note : list() instead of c() avoids numbers being converted to characters
test$group_agg <- as.factor(test$group_agg)

test$prop <- test$nb_fed_day0/test$total
# co-feeding mosquitoes and control are now considered the same group
test$virus_true[test$virus_true == ""] <- "Control"
test$virus_true <- as.factor(test$virus_true)
test <- test[!is.na(test$int_duration),]

# Simple model with duration incorporated as offset (control and concurrent feeding together) ----

# These three formulations are strictly equivalent
# m0 <- glmmTMB(nb_fed_day0/total ~ virus_true*species + offset(log(int_duration)),
#               data = test,
#               weights = total,
#               family = binomial(link = "logit"),
#               REML = F)
# m0 <- glmmTMB(cbind(nb_fed_day0,total - nb_fed_day0) ~ group_agg + offset(log(int_duration)),
#               data = test,
#               family = binomial(link = "logit"),
#               REML = F)
m0 <- glmmTMB(nb_fed_day0/total ~ group_agg + offset(log(int_duration)),
              data = test,
              weights = total,
              family = binomial(link = "logit"),
              REML = F)
MuMIn::AICc(m0) # 315.3186

m1 <- glmmTMB(nb_fed_day0/total ~ group_agg + offset(log(int_duration)),
               data = test,
               family = betabinomial(link = "logit"),
               weights = total,
               dispformula = ~1,
               REML = F)
MuMIn::AICc(m1) # 284.3614
anova(m0,m1) # signif 

m2 <- glmmTMB(nb_fed_day0/total ~ group_agg + offset(log(int_duration)) + (1|ID),
              data = test,
              family = binomial(link = "logit"),
              weights = total,
              REML = F)
MuMIn::AICc(m2) # 292.5736
anova(m0,m2) # signif

m3 <- glmmTMB(nb_fed_day0/total ~ group_agg + offset(log(int_duration)) + (1|ID),
              data = test,
              family = betabinomial(link = "logit"),
              weights = total,
              dispformula = ~1 ,
              REML = F)
MuMIn::AICc(m3) # 286.9377
anova(m2,m3) # signif
anova(m1,m3) # not signif

# saveRDS(m1, "../output/engorgement_models/vector_exposition.rds")
m1 <- readRDS("../output/engorgement_models/vector_exposition.rds")
summary(m1)

comp <- c("Cyno.Zika_virus - Cyno. = 0",
          "Squirrel.Zika_virus - Squirrel. = 0",
          "Squirrel.Dengue_virus - Squirrel. = 0",
          "Cyno.Dengue_virus - Cyno. = 0",
          "Cyno. - Squirrel. = 0",
          "Cyno.Dengue_virus - Squirrel.Dengue_virus = 0",
          "Cyno.Zika_virus - Squirrel.Zika_virus = 0")

tuk = glht(m1, linfct = mcp(group_agg = comp)) 

summary(tuk, test = adjusted("none")) # effect of DENV infection in squirrels + species effect in DENV-infected only
summary(tuk, test = adjusted("fdr")) # same

res <- tidy(confint(tuk))
res$estimate <- exp(res$estimate)
res$conf.low <- exp(res$conf.low)
res$conf.high <- exp(res$conf.high)
# less feeding in infected squirrels than controls + less feeding in DENV-infected squirrels than DENV-infected cynos


# Figure  ----
# Without duration  explicitely represented

newDat <- data.frame(group_agg = unique(test$group_agg))
newDat$int_duration <- median(test$int_duration) # 6 minutes
newDat$total <- 15
# Could not automate this through str_split
newDat$virus_true <- c("Zika virus",
                       "Control",
                       "Zika virus",
                       "Dengue-2 virus",
                       "Dengue-2 virus",
                       "Control")
pred <- predict(m1, newdata = newDat, se.fit = T)
newDat$pred <- plogis(pred$fit)
newDat$lwr <- plogis(pred$fit - 1.96*pred$se.fit)
newDat$upr <- plogis(pred$fit + 1.96*pred$se.fit)


## Model output only (assume equal duration) ----

p0 <- ggplot(newDat, aes(x = as.character(group_agg), y = pred)) +
  geom_errorbar(aes(x = as.character(group_agg),
                    ymin = lwr, ymax = upr,
                    color = virus_true),
                width = 0.1, linewidth = 1.3) +
  geom_point(aes(x = as.character(group_agg),
                 fill = virus_true,
                 shape = virus_true),
             color = "#414341",
             stroke = 1.5,
             size = 6) +
  scale_fill_manual(values = c("Dengue-2 virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue-2 virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#C2c3c9")) +
  scale_shape_manual(values = c("Dengue-2 virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  geom_vline(xintercept = 3.5, alpha = 0.6, linewidth = 1.1, color = "grey") +
  geom_bracket(xmin = "Squirrel.", xmax = "Squirrel.Zika_virus",
               y.position = 1.22, label = "ns", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Cyno.", xmax = "Cyno.Zika_virus",
               y.position = 1.22, label = "ns", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Squirrel.", xmax = "Squirrel.Dengue_virus",
               y.position = 1.08, label = "***", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Cyno.", xmax = "Cyno.Dengue_virus",
               y.position = 1.08, label = "ns", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Cyno.", xmax = "Squirrel.",
               y.position = 1.36, label = "ns", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Cyno.Dengue_virus", xmax = "Squirrel.Dengue_virus",
               y.position = 1.5, label = "***", size = 1.3, label.size = 11) +
  geom_bracket(xmin = "Cyno.Zika_virus", xmax = "Squirrel.Zika_virus",
               y.position = 1.64, label = "ns", size = 1.3, label.size = 11) +
  labs(y = "Predicted probability<br>of engorging",
       x = "",
       fill = "Mosquito\nexposition status",
       color = "Mosquito\nexposition status",
       shape = "Mosquito\nexposition status") +
  ggtitle("A") +
  guides(fill = guide_legend(order = 1,
                             override.aes = list(alpha = 1, size = 5, linewidth = 1)),
         color = guide_legend(order = 1),
         shape = guide_legend(order = 1)) +
  annotate(geom = "text", label = "Cynomolgus\nmacaques",
           x = 2.32, y = 2, size = 13) +
  annotate(geom = "text", label = "Squirrel\nmonkeys",
           x = 5.34, y = 2, size = 13) +
  scale_y_continuous(limits = c(0,2.2),
                     breaks = c(0,0.25,0.5,0.75,1),
                     expand = expansion(add = 0.01)) + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_markdown(size = 38,
                                        margin = margin(r = 30, t = 550),
                                        vjust = 0),
        legend.position = c(1.1,0.27),
        plot.margin = unit(c(0,10,0,0.35), "cm"),
        legend.text = element_text(size = 38),
        legend.title = element_text(size = 38),
        plot.title = element_text(size = 32))

## Raw data ----
p1 <- ggplot(test, aes(x = as.character(group_agg), y = prop)) +
  geom_jitter(aes(fill = virus_true,
                  color = virus_true,
                  shape = virus_true),
              size = 6,
              width = 0.23,
              alpha = 0.7) +
  geom_boxplot(alpha = 0, outlier.shape = NA) +
  scale_fill_manual(values = c("Dengue virus" = "#1c812b",
                               "Zika virus" = "#253dbe",
                               "Control" = "#C2c3c9")) +
  scale_color_manual(values = c("Dengue virus" = "#1c812b",
                                "Zika virus" = "#253dbe",
                                "Control" = "#414341")) +
  scale_shape_manual(values = c("Dengue virus" = 21,
                                "Zika virus" = 22,
                                "Control" = 24)) +
  geom_vline(xintercept = 3.5, alpha = 0.6, linewidth = 1.1, color = "grey") +
  labs(y = "Proportion of<br>engorged mosquitoes",
       x = "",
       fill = "Mosquito status",
       color = "Mosquito status",
       shape = "Mosquito status") +
  guides(fill = guide_legend(order = 1,
                             override.aes = list(alpha = 1, size = 5)),
         color = guide_legend(order = 1),
         shape = guide_legend(order = 1)) +
  scale_y_continuous(limits = c(0,1.05),
                     breaks = c(0,0.25,0.5,0.75,1),
                     expand = expansion(add = 0.01)) + 
  ggtitle("B") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_markdown(size = 38,
                                        margin = margin(r = 30),
                                        vjust = 0),
        legend.position = "none",
        plot.margin = unit(c(0,10,0,0.35), "cm"),
        legend.text = element_text(size = 38),
        legend.title = element_text(size = 38),
        plot.title = element_text(size = 32))

img_cyno <- magick::image_read("../output/figures/outlines/macaque_outline.png")
img_squirrel <- magick::image_read("../output/figures/outlines/squirrel_outline.png")

P <- p0 / p1
P <- P + plot_layout(heights = c(2,1))

p_draw <- ggdraw() + 
  draw_plot(P) +
  draw_image(image = img_cyno, 
             x = 0.15, y = 0.87, scale = 0.1,
             valign = 0, halign = 0) +
  draw_image(image = img_squirrel, 
             x = 0.51, y = 0.87, scale = 0.1,
             valign = 0, halign = 0)

# png(filename = "../output/figures/Figure_2.png",
#     width = 1600, height = 1600)
# print(p_draw)
# dev.off()

pdf(file = "../output/figures/Figure_2.pdf",
    width = 21, height = 21)
print(p_draw)
dev.off()




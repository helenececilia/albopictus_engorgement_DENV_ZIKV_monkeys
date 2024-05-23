## ---------------------------
## Author: Helene Cecilia
##
## Date Created: 2023-02-28

rm(list=ls())

## Loading Packages  ------------------
library(performance) # for check_collinearity
library(ggplot2)
library(tidyverse)
library(BAS)
library(patchwork)
library(cowplot)
library(gghalves)
library(car) # for leveneTest
library(MASS) # for glm.nb and multidimensional scaling (isoMDS, sammon)
library(VGAM) # for vglm
library(lme4) # for lmer
library(lmerTest)
library(modelr)
library(scales) # for trans_breaks
library(MuMIn) # for r.squaredGLMM, AICc
library(lattice) # for qqmath
library(sjPlot) # for plot_model
library(nlme) # for varIdent
library(report)
library(kableExtra)
library(bbmle)
library(clusrank) # wilcoxon test for clustered (non-independent) data
library(DHARMa)
library(effects)
library(glmmTMB)
library(gplots) # for textplot
library(mgcv) # for gam
library(ggtext) # for element_markdown
library(remef) # devtools::install_github("hohenstein/remef")
library(merTools) # for predictInterval
library(gammit) # devtools::install_github("m-clark/gammit") / for predict_gamm
library(wesanderson)
library(gratia) # for confint.gam

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

# Data -----
my_data <- read.csv("../data/Cytokines_and_bites_data.csv")

cyto1 <- my_data[my_data$cytokine == "EGF" & my_data$NHP == "Cyno",]
cyto2 <- my_data[my_data$cytokine == "MIF" & my_data$NHP == "Cyno",]
cyto3 <- my_data[my_data$cytokine == "TGFbeta" & my_data$NHP == "Cyno",]
cyto4 <- my_data[my_data$cytokine == "MCP.1" & my_data$NHP == "Cyno",]

d_cyno <- cyto4
d_cyno$ID <- as.factor(d_cyno$ID)

m_inf <- gam(log10(value) ~ s(bites_of_yesterday, k = 4) + s(cumul_bites_7_previous_days, k = 4) +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")

newDat <- data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                     cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                     ID = d_cyno$ID)
pred_nl <- predict(m_inf, type = "response", newdata = newDat)
# equivalent
# pred_nl2 <- predict_gamm(m_inf, newdata = newDat, re_form = NULL)
newDat$pred_cyto <- pred_nl
newDat$true <- log10(d_cyno$value)
newDat$day <- d_cyno$day

# Data plot ----
bites_only <- read.csv("../data/Control_monkeys_cumulative_bites.csv")
bites_only <- bites_only[bites_only$species == "Cyno",]
bites_no_cyto <- bites_only[bites_only$day %in% c(11,14,21),]
p1 <- ggplot(bites_only) + geom_line(aes(x = day,
                                     y = bites_of_yesterday,
                                     group = ID)) +
  geom_point(aes(x = day,
                y = bites_of_yesterday,
                color = ID),
             alpha = 0.75, size = 6) +
  geom_point(data = bites_no_cyto,
             aes(x = day,
                 y = bites_of_yesterday,
                 fill = ID),
             shape = 21,
             alpha = 0.75, size = 6) +
  facet_wrap(.~ID, nrow = 1) +
  scale_x_continuous(breaks = c(-7,seq(0,7),28),
                     labels = c("-7","0",rep("",6),"7","28")) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme_classic() +
  ggtitle("A") +
  labs(x = "",
       y = "Number of bites<br>in the short-term") +
  theme(legend.position = "none",
        text = element_text(size = 30),
        plot.title = element_text(size = 27),
        axis.title.y = element_markdown())

p2 <- ggplot(bites_only) + geom_line(aes(x = day,
                                     y = cumul_bites_7_previous_days,
                                     group = ID)) +
  geom_point(aes(x = day,
                y = cumul_bites_7_previous_days,
                color = ID),
             alpha = 0.75, size = 6) +
  geom_point(data = bites_no_cyto,
             aes(x = day,
                 y = cumul_bites_7_previous_days,
                 fill = ID),
             shape = 21,
             alpha = 0.75, size = 6) +
  scale_x_continuous(breaks = c(-7,seq(0,7),28),
                     labels = c("-7","0",rep("",6),"7","28")) +
  facet_wrap(.~ID, nrow = 1) +
  ggtitle("B") +
  theme_classic() +
  labs(x = "Day post infection",
       y = "Number of bites<br>in the long-term") +
  theme(legend.position = "none",
        text = element_text(size = 30),
        plot.title = element_text(size = 27),
        axis.title.y = element_markdown())

png(filename = "../output/figures/Figure_S3.png",
    width = 1266, height = 800)
plot(p1 / p2)
dev.off()

# Model fit plot ----
p3 <- ggplot(newDat) +
  geom_line(aes(x = day,y = pred_cyto,
                color = ID,
                group = ID)) +
  geom_point(aes(x = day,
                 y = true,
                 color = ID),
             alpha = 0.75, size = 6) +
  facet_wrap(.~ID, ncol = 1) +
  # ggtitle("A") +
  scale_x_continuous(breaks = c(-7,seq(0,7),28),
                     labels = c("-7","0",rep("",6),"7","28")) +
  theme_classic() +
  labs(x = "Day post infection",
       y = "Cytokine concentration (log<sub>10</sub> pg/\u03BCl)") +
  theme(legend.position = "none",
        text = element_text(size = 30),
        plot.title = element_text(size = 27),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())

p4 <- ggplot() +
  geom_line(data = newDat,
            aes(x = bites_of_yesterday, y = pred_cyto,
                group = ID, color = ID)) +
  geom_point(data = d_cyno,
             aes(x = bites_of_yesterday, y = log10(value),
                 color = ID),
             alpha = 0.75, size = 6) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  facet_wrap(.~ID, ncol = 1) +
  # ggtitle("B") +
  theme_classic() +
  labs(x = "Number of bites<br>in the short-term",
       y = "")+
  theme(legend.position = "none",
        text = element_text(size = 30),
        plot.title = element_text(size = 27),
        axis.title.x = element_markdown())

p5 <- ggplot() +
  geom_line(data = newDat,
            aes(x = cumul_bites_7_previous_days, y = pred_cyto,
                group = ID, color = ID)) +
  geom_point(data = d_cyno,
             aes(x = cumul_bites_7_previous_days, y = log10(value),
                 color = ID),
             alpha = 0.75, size = 6) +
  facet_wrap(.~ID, ncol = 1) +
  # ggtitle("C") +
  theme_classic() +
  labs(x = "Number of bites<br>in the long-term",
       y = "")+
  theme(legend.position = "none",
        text = element_text(size = 30),
        plot.title = element_text(size = 27),
        axis.title.x = element_markdown())

p <- (p3 | p4 | p5) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 25))

png(filename = "../output/figures/Figure_S4.png",
    width = 1900, height = 1000) # cyto1
# png(filename = "../output/figures/Figure_S5.png",
#     width = 1900, height = 1000) # cyto2
plot(p)
dev.off()

# for supplementary cyto (effect long term only)
p <- (p3 | p5) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 25))

png(filename = "../output/figures/Figure_S7.png",
    width = 1266, height = 1000) # cyto3
# png(filename = "../output/figures/Figure_S8.png",
#     width = 1266, height = 1000) # cyto 4
plot(p)
dev.off()


# glmmTMB version ----
pdf(file = "../output/cytokines_analysis/report_repeated_bites_controls_glmmTMB_no_LOD_AICc.pdf")
for(c in unique(data$cytokine)){
  df <- my_data[my_data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Bites in cyno", cex = 2)
  d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
  nb_excl <- length(df[df$value == df$LOD & df$NHP == "Cyno",]$value)
  # browser()
  textplot(paste0("Nb excluded (LOD): ", nb_excl, "\nNb remaining: ", length(d_cyno$value)),
           cex = 2, mar = c(5,5,5,5))
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ bites_of_yesterday + cumul_bites_7_previous_days +
                       (1|ID),
                     REML = FALSE,
                     data = d_cyno)
    simulateResiduals(m_inf, plot = T)
    plot(allEffects(m_inf, partial.residuals = T))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                                                                      cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                                                                      ID = d_cyno$ID,
                                                                      day = d_cyno$day))
    pred <- data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                       cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = bites_of_yesterday, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = bites_of_yesterday, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites yesterday",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Bites in squirrel", cex = 2)
  d_sq <- df[df$NHP == "Squirrel",]
  nb_excl <- length(d_sq[d_sq$value == d_sq$LOD,]$value)
  d_sq <- d_sq[d_sq$value != d_sq$LOD,]
  textplot(paste0("Nb excluded (LOD): ", nb_excl, "\nNb remaining: ", length(d_sq$value)),
           cex = 2, mar = c(5,5,5,5))  
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ cumul_bites_2_previous_days + cumul_bites_7_previous_days +
                       (1|ID),
                     REML = FALSE,
                     data = d_sq)
    simulateResiduals(m_inf, plot = T)
    plot(allEffects(m_inf, partial.residuals = T))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                                                                      cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                                                                      ID = d_sq$ID,
                                                                      day = d_sq$day))
    pred <- data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                       cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_2_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_2_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 2 previous days",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
}
dev.off()


pdf(file = "../output/cytokines_analysis/report_repeated_bites_controls_glmmTMB_with_LOD_AICc.pdf")
for(c in unique(data$cytokine)){
  df <- my_data[my_data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Bites in cyno", cex = 2)
  d_cyno <- df[df$NHP == "Cyno",]
  # browser()
  textplot(paste0("Nb obs: ", length(d_cyno$value)),
           cex = 2, mar = c(5,5,5,5))
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ bites_of_yesterday + cumul_bites_7_previous_days +
                       (1|ID),
                     REML = FALSE,
                     data = d_cyno)
    simulateResiduals(m_inf, plot = T)
    plot(allEffects(m_inf, partial.residuals = T))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                                                                      cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                                                                      ID = d_cyno$ID,
                                                                      day = d_cyno$day))
    pred <- data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                       cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = bites_of_yesterday, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = bites_of_yesterday, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites yesterday",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Bites in squirrel", cex = 2)
  d_sq <- df[df$NHP == "Squirrel",]
  
  textplot(paste0("Nb obs: ", length(d_sq$value)),
           cex = 2, mar = c(5,5,5,5))  
  tryCatch({
    m_inf <- glmmTMB(log10(value) ~ cumul_bites_2_previous_days + cumul_bites_7_previous_days +
                       (1|ID),
                     REML = FALSE,
                     data = d_sq)
    simulateResiduals(m_inf, plot = T)
    plot(allEffects(m_inf, partial.residuals = T))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                                                                      cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                                                                      ID = d_sq$ID,
                                                                      day = d_sq$day))
    pred <- data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                       cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_2_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_2_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 2 previous days",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
}
dev.off()


# gam version ----
pdf(file = "../output/cytokines_analysis/report_repeated_bites_controls_GAM_no_LOD_AICc.pdf")
for(c in unique(data$cytokine)){
  df <- my_data[my_data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Bites in cyno", cex = 2)
  d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
  nb_excl <- length(df[df$value == df$LOD & df$NHP == "Cyno",]$value)
  # browser()
  textplot(paste0("Nb excluded (LOD): ", nb_excl, "\nNb remaining: ", length(d_cyno$value)),
           cex = 2, mar = c(5,5,5,5))
  d_cyno$ID <- as.factor(d_cyno$ID)
  tryCatch({
    m_inf <- gam(log10(value) ~ s(bites_of_yesterday, k = 4) + s(cumul_bites_7_previous_days, k = 4) +
                   s(ID, bs = "re", k = 2),
                 data = d_cyno,
                 method = "ML")
    # simulateResiduals(m_inf, plot = T)
    plot(m_inf, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
         col = c("black", rep("red", length(m_inf$residuals))))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                                                                      cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                                                                      ID = d_cyno$ID))
    pred <- data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                       cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = bites_of_yesterday, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = bites_of_yesterday, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites yesterday",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(gam.check(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Bites in squirrel", cex = 2)
  d_sq <- df[df$NHP == "Squirrel",]
  nb_excl <- length(d_sq[d_sq$value == d_sq$LOD,]$value)
  d_sq <- d_sq[d_sq$value != d_sq$LOD,]
  d_sq$ID <- as.factor(d_sq$ID)
  textplot(paste0("Nb excluded (LOD): ", nb_excl, "\nNb remaining: ", length(d_sq$value)),
           cex = 2, mar = c(5,5,5,5))  
  tryCatch({
    # Because of the sampling scheme of squirrels, bites_of_yesterday is useless
    m_inf <- gam(log10(value) ~ s(cumul_bites_2_previous_days, k = 4) + s(cumul_bites_7_previous_days, k = 4) +
                   s(ID, bs = "re", k = 2),
                 data = d_sq,
                 method = "ML")
    # simulateResiduals(m_inf, plot = T)
    plot(m_inf, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
         col = c("black", rep("red", length(m_inf$residuals))))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                                                                      cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                                                                      ID = d_sq$ID))
    pred <- data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                       cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_2_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_2_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 2 previous days",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(gam.check(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
}
dev.off()

pdf(file = "../output/cytokines_analysis/report_repeated_bites_controls_GAM_with_LOD_AICc.pdf")
for(c in unique(data$cytokine)){
  df <- my_data[my_data$cytokine == c,]
  textplot(c, cex = 2)
  #########
  textplot("Bites in cyno", cex = 2)
  d_cyno <- df[df$NHP == "Cyno",]
  d_cyno$ID <- as.factor(d_cyno$ID)
  textplot(paste0("Nb obs: ", length(d_cyno$value)),
           cex = 2, mar = c(5,5,5,5))
  tryCatch({
    m_inf <- gam(log10(value) ~ s(bites_of_yesterday, k = 4) + s(cumul_bites_7_previous_days, k = 4) +
                   s(ID, bs = "re", k = 2),
                 data = d_cyno,
                 method = "ML")
    # simulateResiduals(m_inf, plot = T)
    plot(m_inf, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
         col = c("black", rep("red", length(m_inf$residuals))))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                                                                      cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                                                                      ID = d_cyno$ID))
    pred <- data.frame(bites_of_yesterday = d_cyno$bites_of_yesterday,
                       cumul_bites_7_previous_days = d_cyno$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = bites_of_yesterday, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = bites_of_yesterday, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites yesterday",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_cyno,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(gam.check(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
  
  textplot("Bites in squirrel", cex = 2)
  d_sq <- df[df$NHP == "Squirrel",]
  
  d_sq$ID <- as.factor(d_sq$ID)
  textplot(paste0("Nb obs: ", length(d_sq$value)),
           cex = 2, mar = c(5,5,5,5))  
  tryCatch({
    # Because of the sampling scheme of squirrels, bites_of_yesterday is useless
    m_inf <- gam(log10(value) ~ s(cumul_bites_2_previous_days, k = 4) + s(cumul_bites_7_previous_days, k = 4) +
                   s(ID, bs = "re", k = 2),
                 data = d_sq,
                 method = "ML")
    # simulateResiduals(m_inf, plot = T)
    plot(m_inf, pages = 0, residuals = T, pch = 20, lwd = 1.8, cex = 0.7,
         col = c("black", rep("red", length(m_inf$residuals))))
    
    pred_nl <- predict(m_inf, type = "response", newdata = data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                                                                      cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                                                                      ID = d_sq$ID))
    pred <- data.frame(cumul_bites_2_previous_days = d_sq$cumul_bites_2_previous_days,
                       cumul_bites_7_previous_days = d_sq$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl)

    p1 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_2_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_2_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 2 previous days",
           y = "Cytokine concentration (log10)")

    p2 <- ggplot() +
      geom_line(data = pred,
                aes(x = cumul_bites_7_previous_days, y = pred_cyto)) +
      geom_point(data = d_sq,
                 aes(x = cumul_bites_7_previous_days, y = log10(value)),
                 alpha = 0.75, size = 3, color = "orange") +
      theme_classic() +
      labs(x = "Nb bites 7 previous days",
           y = "Cytokine concentration (log10)")
    plot(p1|p2)
    
    textplot(capture.output(gam.check(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(check_collinearity(m_inf)), mar = c(5,5,5,5))
    textplot(capture.output(summary(m_inf)), mar = c(5,5,5,5))
    textplot(paste0("AICc ",capture.output(MuMIn::AICc(m_inf))), mar = c(5,5,5,5))
  }, error=function(e){textplot(capture.output(cat(c, "ERROR :",conditionMessage(e), "\n")), mar = c(5,5,5,5))})
}
dev.off()



# Plots ----
# Cyno : short term effect : EGF and MIF gam
# Cyno : long term effect : EGF and MIF gam
# Cyno : TGF beta & MCP.1 : long term 

# Gam models with significant effects ----
pal <- c("NV259" = "#c52415",
         "NV289" = "#7112ac",
         "UG171" = "#7b5c0d",
         "UG253A" = "#f6b81a")

## Short term ----
c <- "EGF" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)

m_inf <- gam(log10(value) ~ bites_of_yesterday +
               cumul_bites_7_previous_days +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   log10(value) ~ bites_of_yesterday + cumul_bites_7_previous_days + 
#   s(ID, bs = "re", k = 2)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  1.994908   0.049363  40.413  < 2e-16 ***
#   bites_of_yesterday           0.019460   0.004632   4.201 0.000215 ***
#   cumul_bites_7_previous_days -0.006161   0.001149  -5.363 8.03e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(ID) 2.596      3 8.062 0.000166 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.619   Deviance explained = 66.9%
# -ML = -31.996  Scale est. = 0.0084674  n = 36

min_ <- min(d_cyno$bites_of_yesterday)
max_ <- max(d_cyno$bites_of_yesterday)
newDat <- data.frame(bites_of_yesterday = seq(min_,max_,0.1))
newDat$cumul_bites_7_previous_days <- 10 
newDat$ID <- "NV259" # needed to run but not used (re_form = NA)
pred_nl_den <- predict_gamm(m_inf, type = "response",
                       newdata = newDat, se = T,
                       re_form = NA)
pred_den <- data.frame(bites = newDat$bites_of_yesterday,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$bites_of_yesterday)
  max_ <- max(df$bites_of_yesterday)
  new <- data.frame(bites_of_yesterday = seq(min_,max_,0.1))
  new$cumul_bites_7_previous_days <- 10 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$bites_of_yesterday,
                       pred_cyto = pred_nl_den$prediction,
                       ID = newDat$ID)

p_egf <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                              linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  scale_color_manual(values = pal) +
  labs(y = "Concentration (log<sub>10</sub> pg/\u03BCl)",
       x = "Number of bites in the short-term",
       color = "") +
  ggtitle(paste0("A - ",c)) +
  scale_y_continuous(limits = c(1.5,2.35),
                     breaks = c(1.5,1.7,1.9,2.1,2.3),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10),
                     expand = expansion(add = 0.12)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = c(0.8,0.33),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28),
        plot.margin = margin(r = 10))

c <- "MIF" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)
m_inf <- gam(log10(value) ~ s(bites_of_yesterday, k = 4) +
               cumul_bites_7_previous_days +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   log10(value) ~ s(bites_of_yesterday, k = 4) + cumul_bites_7_previous_days + 
#   s(ID, bs = "re", k = 2)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  2.474104   0.075431  32.800  < 2e-16 ***
#   cumul_bites_7_previous_days -0.009529   0.002281  -4.177 0.000234 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(bites_of_yesterday) 1.736  2.052 6.168 0.00477 **
#   s(ID)                 2.202  3.000 3.533 0.00718 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.526   Deviance explained = 59.3%
# -ML = -6.7505  Scale est. = 0.033338  n = 36
min_ <- min(d_cyno$bites_of_yesterday)
max_ <- max(d_cyno$bites_of_yesterday)
newDat <- data.frame(bites_of_yesterday = seq(min_,max_,0.1))
newDat$cumul_bites_7_previous_days <- 10 
newDat$ID <- "NV259" 
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA)
pred_den <- data.frame(bites = newDat$bites_of_yesterday,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$bites_of_yesterday)
  max_ <- max(df$bites_of_yesterday)
  new <- data.frame(bites_of_yesterday = seq(min_,max_,0.1))
  new$cumul_bites_7_previous_days <- 10 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$bites_of_yesterday,
                        pred_cyto = pred_nl_den$prediction,
                        ID = newDat$ID)

p_mif <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                              linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  labs(y = "", 
       x = "Number of bites in the short-term",
       color = "") +
  ggtitle(paste0("B - ",c)) +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(1.57,2.9),
                     breaks = c(1.7,2.0,2.3,2.6,2.9),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10),
                     expand = expansion(add = 0.1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = "none",
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))

p_short <- p_egf | p_mif

## Long term ----
c <- "EGF" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)

m_inf <- gam(log10(value) ~ bites_of_yesterday +
               cumul_bites_7_previous_days+
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)
min_ <- min(d_cyno$cumul_bites_7_previous_days)
max_ <- max(d_cyno$cumul_bites_7_previous_days)
newDat <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
newDat$bites_of_yesterday <- 0 
newDat$ID <- "NV259" 
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA)
pred_den <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$cumul_bites_7_previous_days)
  max_ <- max(df$cumul_bites_7_previous_days)
  new <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
  new$bites_of_yesterday <- 0 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                        pred_cyto = pred_nl_den$prediction,
                        ID = newDat$ID)

p_egf <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                              linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  labs(y = "Concentration (log<sub>10</sub> pg/\u03BCl)",
       x = "Number of bites in the long-term",
       color = "") +
  ggtitle(paste0("C - ",c)) +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(1.5,2.35),
                     breaks = c(1.5,1.7,1.9,2.1,2.3),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(limits = c(0,46),
                     breaks = c(0,5,10,15,20,25,30,35,40,45),
                     expand = expansion(add = 0.12)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = "none", #c(0.33,0.25),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28),
        plot.margin = margin(r = 10))

c <- "MIF" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)
m_inf <- gam(log10(value) ~ s(bites_of_yesterday, k = 4) +
               s(cumul_bites_7_previous_days, k = 4) +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)
min_ <- min(d_cyno$cumul_bites_7_previous_days)
max_ <- max(d_cyno$cumul_bites_7_previous_days)
newDat <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
newDat$bites_of_yesterday <- 0 
newDat$ID <- "NV259" 
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA)
pred_den <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$cumul_bites_7_previous_days)
  max_ <- max(df$cumul_bites_7_previous_days)
  new <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
  new$bites_of_yesterday <- 0 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                        pred_cyto = pred_nl_den$prediction,
                        ID = newDat$ID)

p_mif <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                              linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  labs(y = "", 
       x = "Number of bites in the long-term",
       color = "") +
  ggtitle(paste0("D - ",c)) +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(1.57,2.9),
                     breaks = c(1.7,2.0,2.3,2.6,2.9),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(limits = c(0,46),
                     breaks = c(0,5,10,15,20,25,30,35,40,45),
                     expand = expansion(add = 0.12)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = "none",
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))

p_long <- p_egf | p_mif


p <- p_short / p_long

p_draw <- ggdraw(p) + 
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = -0.41, y = 0.4, scale = 0.07) +
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = -0.41, y = -0.1, scale = 0.07) +
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = 0.08, y = 0.4, scale = 0.07) +
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = 0.08, y = -0.1, scale = 0.07) 

# png(filename = "../output/figures/Figure_6.png", width = 1600, height = 1000)
# print(p_draw)
# dev.off()

cairo_pdf(filename = "../output/figures/Figure_6.pdf",
          width = 21, height = 13)
print(p_draw)
dev.off()

# Supplementary ----
c <- "TGFbeta" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)

m_inf <- gam(log10(value) ~ bites_of_yesterday +
               cumul_bites_7_previous_days +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   log10(value) ~ bites_of_yesterday + cumul_bites_7_previous_days + 
#   s(ID, bs = "re", k = 2)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  1.927755   0.044615  43.208  < 2e-16 ***
#   bites_of_yesterday           0.006916   0.004890   1.414  0.16739    
# cumul_bites_7_previous_days -0.003627   0.001212  -2.992  0.00543 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value   
# s(ID) 2.361      3 4.58 0.00282 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.403   Deviance explained = 47.8%
# -ML = -30.99  Scale est. = 0.0094418  n = 36


# long term

min_ <- min(d_cyno$cumul_bites_7_previous_days)
max_ <- max(d_cyno$cumul_bites_7_previous_days)
newDat <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
newDat$bites_of_yesterday <- 0 
newDat$ID <- "NV259" 
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA)
pred_den <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$cumul_bites_7_previous_days)
  max_ <- max(df$cumul_bites_7_previous_days)
  new <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
  new$bites_of_yesterday <- 0 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                        pred_cyto = pred_nl_den$prediction,
                        ID = newDat$ID)

p_tgf_long <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                                   linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  labs(y = "Concentration (log<sub>10</sub> pg/\u03BCl)",
       x = "Number of bites in the long-term",
       color = "") +
  scale_color_manual(values = pal) +
  ggtitle(paste0("A - ",c)) +
  scale_y_continuous(limits = c(1.6,2.1),
                     breaks = c(1.6,1.7,1.8,1.9,2.0,2.1),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(limits = c(0,46),
                     breaks = c(0,5,10,15,20,25,30,35,40,45),
                     expand = expansion(add = 0.12)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = c(0.75,0.8),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))

c <- "MCP.1" 
df <- my_data[my_data$cytokine == c,]
d_cyno <- df[df$value != df$LOD & df$NHP == "Cyno",]
d_cyno$ID <- as.factor(d_cyno$ID)

m_inf <- gam(log10(value) ~ bites_of_yesterday +
               cumul_bites_7_previous_days +
               s(ID, bs = "re", k = 2),
             data = d_cyno,
             method = "ML")
MuMIn::AICc(m_inf)
summary(m_inf)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   log10(value) ~ bites_of_yesterday + cumul_bites_7_previous_days + 
#   s(ID, bs = "re", k = 2)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  2.743543   0.125409  21.877  < 2e-16 ***
#   bites_of_yesterday           0.007956   0.005265   1.511  0.14121    
# cumul_bites_7_previous_days -0.003885   0.001307  -2.972  0.00577 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(ID) 2.938      3 58.49  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.837   Deviance explained =   86%
# -ML = -23.593  Scale est. = 0.010935  n = 36


# long term

min_ <- min(d_cyno$cumul_bites_7_previous_days)
max_ <- max(d_cyno$cumul_bites_7_previous_days)
newDat <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
newDat$bites_of_yesterday <- 0 
newDat$ID <- "NV259" 
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = T,
                            re_form = NA)
pred_den <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                       pred_cyto = pred_nl_den$prediction,
                       lwr = pred_nl_den$prediction + qnorm(0.025)*pred_nl_den$se,
                       upr = pred_nl_den$prediction + qnorm(0.975)*pred_nl_den$se)

newDat <- data.frame()
for(i in unique(d_cyno$ID)){
  df <- d_cyno[d_cyno$ID == i,]
  min_ <- min(df$cumul_bites_7_previous_days)
  max_ <- max(df$cumul_bites_7_previous_days)
  new <- data.frame(cumul_bites_7_previous_days = seq(min_,max_,0.1))
  new$bites_of_yesterday <- 0 
  new$ID <- i
  newDat <- rbind(newDat,new)
}
pred_nl_den <- predict_gamm(m_inf, type = "response",
                            newdata = newDat, se = F,
                            re_form = NULL) # NULL = with ranef / NA = no ranef
pred_den2 <- data.frame(bites = newDat$cumul_bites_7_previous_days,
                        pred_cyto = pred_nl_den$prediction,
                        ID = newDat$ID)

p_mcp_long <- ggplot() + geom_line(data = pred_den, aes(x = bites, y = pred_cyto),
                                   linewidth = 1.5) +
  geom_line(data = pred_den2, aes(x = bites, y = pred_cyto, group = ID, col = ID),
            linewidth = 1.1, alpha = 0.75) +
  geom_ribbon(data = pred_den, aes(x = bites, ymin = lwr, ymax = upr),
              fill = "lightgrey", alpha = 0.3) +
  labs(y = "", 
       x = "Number of bites in the long-term",
       color = "") +
  scale_color_manual(values = pal) +
  ggtitle(paste0("B - ",c)) +
  scale_y_continuous(limits = c(2.2,3.21),
                     breaks = c(2.2,2.4,2.6,2.8,3.0,3.2),
                     expand = expansion(add = 0.01)) +
  scale_x_continuous(limits = c(0,46),
                     breaks = c(0,5,10,15,20,25,30,35,40,45),
                     expand = expansion(add = 0.12)) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title.y = element_markdown(size = 30,
                                        margin = margin(r = 10)),
        axis.title.x = element_text(size = 30),
        legend.position = "none", #c(0.33,0.25),
        legend.text = element_text(size = 28),
        plot.title = element_text(size = 28))

p_suppl <- p_tgf_long | p_mcp_long

p_draw <- ggdraw(p_suppl) + 
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = -0.41, y = -0.25, scale = 0.12) +
  draw_image(image = "../output/figures/outlines/macaque_outline.png",
             x = 0.08, y = -0.25, scale = 0.12) 

png(filename = "../output/figures/Figure_S3.png", width = 1600, height = 500)
plot(p_draw)
dev.off()

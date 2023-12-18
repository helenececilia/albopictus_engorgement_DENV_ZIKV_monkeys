## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Helene Cecilia
##
## Date Created: 2023-06-27

rm(list=ls())

## Loading Packages  ------------------
library(janitor)
library(dplyr)
library(ggplot2)
library(patchwork)

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

df <- read.csv("../data/Approach_rate_060623.csv", sep = "\t", dec = ".")
df <- df %>% clean_names()

stat <- df %>% dplyr::group_by(collection_window) %>%
  dplyr::summarise(prop_non_zero = sum(females_ae_albopictus != 0)/n(),
                   mean_tot = mean(females_ae_albopictus),
                   mean_sub = mean(females_ae_albopictus[females_ae_albopictus != 0]),
                   max = max(females_ae_albopictus)) %>%
  ungroup()

# 12:00-13:00 and 15:00-16:00 are close
# 12-13 has a higher mean in non-zero values
# 15-16 has a higher proportion of non-zero values (and therefore a higher mean of all values)
# both have the same max

stat <- df %>% dplyr::group_by(ndbi) %>%
  dplyr::summarise(prop_non_zero = sum(females_ae_albopictus != 0)/n(),
                   mean_tot = mean(females_ae_albopictus),
                   mean_sub = mean(females_ae_albopictus[females_ae_albopictus != 0]),
                   max = max(females_ae_albopictus)) %>%
  ungroup()

# same conundrum for ndbi 1 and 2
# 1 has more non-zero values (and therefore a higher mean of all values)
# 2 has a higher mean in non-zero values
# both have the same max

stat <- df %>% dplyr::group_by(ndbi,collection_window) %>%
  dplyr::summarise(prop_non_zero = sum(females_ae_albopictus != 0)/n(),
                   mean_tot = mean(females_ae_albopictus),
                   mean_sub = mean(females_ae_albopictus[females_ae_albopictus != 0]),
                   max = max(females_ae_albopictus),
                   nb = n()) %>%
  ungroup()

# 1 x 12-13 seems to be optimal

# coeff of variation per collector for this particular setting
# small number of collections per collector, meaningless
stat <- df[df$ndbi == 1 & df$collection_window == "12:00-13:00",] %>% dplyr::group_by(collector) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus)) %>%
  ungroup()

# restrain to ndbi 1 and 2, 12:13 or 15:16
optimum <- df[df$ndbi %in% c(1,2) & df$collection_window %in% c("12:00-13:00","15:00-16:00"),]
stat <- df[df$ndbi %in% c(1,2) & df$collection_window %in% c("12:00-13:00","15:00-16:00"),] %>%
  dplyr::group_by(collector) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus),
                   nb = n()) %>%
  ungroup()
mean(c(0,2,2,1.24,1.92)) # assigning 0 to the collector with always 0 mosquitoes in several attempts
# collector     cv    nb
# <chr>      <dbl> <int>
# AGDS       NA        1
# AH          2        4
# CR          2        4
# EH          1.24    10
# IP        NaN        6
# TN          1.92     7

optimum$ID <- ""
optimum$ID[optimum$collector == "IP"] <- "A"
optimum$ID[optimum$collector == "AGDS"] <- "B"
optimum$ID[optimum$collector == "AH"] <- "C"
optimum$ID[optimum$collector == "CR"] <- "D"
optimum$ID[optimum$collector == "EH"] <- "E"
optimum$ID[optimum$collector == "TN"] <- "F"


p <- ggplot(optimum) +
  geom_boxplot(aes(x = ID, y = females_ae_albopictus,
                  group = ID),
               color = "darkgrey",
               fill = NA,
               outlier.shape = NA) +
  geom_jitter(aes(x = ID, y = females_ae_albopictus,
                  shape = ID),
              height = 0.05, width = 0.2,
              alpha = 0.7,
              size = 5) +
  labs(shape = "Collector",
       x = "") +
  scale_y_continuous(name = "Number approached female Ae. albopictus",
                     limits = c(-0.05,6.2)) + 
  scale_shape_manual(values = c("A" = 0,
                                "B" = 16,
                                "C" = 2,
                                "D" = 18,
                                "E" = 8,
                                "F" = 15)) +
  annotate(geom = "text", label = "B",
           x = 0.5, y = 6.15, size = 11) + # instead of title or tag because of font
  # theme_cowplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 26),
        axis.title.y = element_text(size = 28,
                                    margin = margin(r = 10)),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        legend.position = "right")

img_human <- magick::image_read("~/Documents/POSTDOC/Presentations/Images/homme.png")

p_draw <- ggdraw() + 
  draw_plot(p) +
  draw_image(image = img_human, 
             x = 0.1, y = 0.8, scale = 0.15,
             valign = 0, halign = 0)

png(filename = "../output/figures/approach_optimum_conditions.png",
    height = 800, width = 1600)
plot(p_draw)
dev.off()

# ndbi_labs <- c("NDBI 1","NDBI 2")
# names(ndbi_labs) <- c("1","2")
# 
# time_labs <- c("From 12:00 to 13:00","From 15:00 to 16:00")
# names(time_labs) <- c("12:00-13:00","15:00-16:00")
#
# p <- ggplot(optimum) +
#   geom_jitter(aes(x = ID, y = females_ae_albopictus,
#                  shape = ID),
#              height = 0.1, width = 0.45,
#              alpha = 0.7,
#              size = 3) +
#   facet_grid(ndbi ~ collection_window,
#              labeller = labeller(ndbi = ndbi_labs,
#                                  collection_window = time_labs)) +
#   labs(shape = "Collector", y = "Nb approached females Ae. albopictus",
#        x = "") +
#   scale_shape_manual(values = c("#A" = 0,
#                                 "#B" = 16,
#                                 "#C" = 2,
#                                 "#D" = 18,
#                                 "#E" = 8,
#                                 "#F" = 15)) +
#   ggtitle("B") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 23),
#         axis.title = element_text(size = 23),
#         strip.text = element_text(size = 22),
#         plot.title = element_text(size = 23),
#         legend.title = element_text(size = 23),
#         legend.text = element_text(size = 22),
#         legend.position = "right")
# ndbi 2 does not have much mosq collected for these collection windows...



# ----
stat <- df %>% dplyr::group_by(collector,ndbi) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus)) %>%
  ungroup()
# # A tibble: 16 × 3
# collector  ndbi     cv
# <chr>     <int>  <dbl>
# AGDS          1  NA   
# AH            1   1.40
# AH            2   2.83
# AH            3 NaN   
# CR            1   1.10
# CR            2 NaN   
# CR            3 NaN   
# EH            1   1.33
# EH            2   3.55
# EH            3 NaN   
# IP            1   1.85
# IP            2 NaN   
# IP            3 NaN   
# TN            1 NaN   
# TN            2   1.35
# TN            3 NaN

stat <- df %>% dplyr::group_by(collector) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus)) %>%
  ungroup()
# # A tibble: 6 × 2
# collector    cv
# <chr>     <dbl>
# AGDS      NA   
# AH         2.62
# CR         2.30
# EH         2.46
# IP         3.39
# TN         2.81

stat <- df %>% dplyr::group_by(ndbi) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus)) %>%
  ungroup()
# # A tibble: 3 × 2
# ndbi     cv
# <int>  <dbl>
#   1   1.61
#   2   3.22
#   3 NaN  
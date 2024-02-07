## ---------------------------
## Author: Helene Cecilia
##
## Date Created: 2023-06-27

rm(list=ls())

## Loading Packages  ------------------
library(janitor)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot) # for draw_image
library(magick) # for draw_image 

## Set Work Directory ------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set to source file location
getwd()

## Load source files ------------------
#source('.R')

## Global command
`%notin%` <- Negate(`%in%`)

## -------------------------------

df <- read.csv("../data/Approach_rate_Hernandez_Acosta.csv", sep = ",", dec = ".")
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

# ndbi 1 and 2 are close
# 1 has more non-zero values (and therefore a higher mean of all values)
# 2 has a higher mean in non-zero values
# both have the same max

# restrain to ndbi 1 and 2, 12:13 or 15:16
optimum <- df[df$ndbi %in% c(1,2) & df$collection_window %in% c("12:00-13:00","15:00-16:00"),]
stat <- optimum %>%
  dplyr::group_by(collector) %>%
  dplyr::summarise(cv = sd(females_ae_albopictus)/mean(females_ae_albopictus),
                   nb = n()) %>%
  ungroup()

mean(c(0,2,2,1.24,1.92)) # assigning 0 to the collector with always 0 mosquitoes in several attempts
# collector     cv    nb
 # A         NaN        6
 # B          NA        1
 # C           2        4
 # D           2        4
 # E           1.24    10
 # F           1.92     7


p <- ggplot(optimum) +
  geom_boxplot(aes(x = id, y = females_ae_albopictus,
                  group = id),
               color = "darkgrey",
               fill = NA,
               outlier.shape = NA) +
  geom_jitter(aes(x = id, y = females_ae_albopictus,
                  shape = id),
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

img_human <- magick::image_read("../output/figures/outlines/human_outline.png")

p_draw <- ggdraw() + 
  draw_plot(p) +
  draw_image(image = img_human, 
             x = 0.1, y = 0.8, scale = 0.15,
             valign = 0, halign = 0)

png(filename = "../output/figures/Figure_S2_B.png",
    height = 800, width = 1600)
plot(p_draw)
dev.off()

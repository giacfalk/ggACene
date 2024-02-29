
#######

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(haven)
library(readstata13)
library(stringr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(ResourceSelection)
library(multiwayvcov)
library(msm) # https://stats.oarc.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
library(margins)
library(texreg)
library(xtable)
library(stargazer)
library(effects)
library(splines)
library(sjPlot)
library(MASS)
library(ISLR)
library(car)
library(ROCR)
library(lspline)
library(sf)
library(pROC)
library(pbapply)

setwd(wd)

# descriptive statistics
#######

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

#######

shape_ac <- shape_ac %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))
shape_ely_diff <- shape_ely_diff %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

########

remotes::install_github("chris-prener/biscale")

library(biscale)

shape_ac_v <- shape_ac
shape_ac_v$var1 <-exp(shape_ac_v$GDP_SSP2_2010)
shape_ac_v$var2 <- shape_ac_v$`CDD_245_2020_GFDL-ESM4`

shape_ac_v <- shape_ac_v %>% group_by(region) %>% mutate(var1=ntile(var1, 5))

shape_ac_v <- dplyr::select(shape_ac_v, var1, var2)
shape_ac_v <- na.omit(shape_ac_v)

data <- bi_class(shape_ac_v, x = var1, y = var2, style = "quantile", dim = 5)

###

map <- ggplot() +
  geom_sf(data = shape_ac, mapping = aes(fill = region), color = "transparent", size = 0, show.legend = T) +
  bi_theme()

ggsave("regions_map.png", map, height=3, width = 10, scale=1.5)

###

custom_pal <- c(
  "1-1" = "#d3d3d3", # low x, low y
  "2-1" = "#b6cdcd",
  "3-1" = "#97c5c5",
  "4-1" = "#75bebe",
  "5-1" = "#52b6b6", # high x, low y
  "1-2" = "#cab6c5",
  "2-2" = "#aeb0bf",
  "3-2" = "#91aab9",
  "4-2" = "#70a4b2",
  "5-2" = "#4e9daa",
  "1-3" = "#c098b9",
  "2-3" = "#a593b3",
  "3-3" = "#898ead",
  "4-3" = "#6b89a6",
  "5-3" = "#4a839f",
  "1-4" = "#b77aab",
  "2-4" = "#9e76a6",
  "3-4" = "#8372a0",
  "4-4" = "#666e9a",
  "5-4" = "#476993",
  "1-5" = "#ad5b9c", # low x, high y
  "2-5" = "#955898",
  "3-5" = "#7c5592",
  "4-5" = "#60528d",
  "5-5" = "#434e87" # high x, high y
)

map <- ggplot() +
  geom_sf(data = data, mapping = aes(fill = bi_class), color = "transparent", size = 0, show.legend = FALSE) +
  bi_scale_fill(pal = custom_pal, dim = 5) +
  labs(
    title = "CDDs and Income distribution",
    caption = "Income quintiles are calculated at the macroregion level."
  ) +
  bi_theme()

legend <- bi_legend(pal = custom_pal,
                    dim = 5,
                    xlab = "Higher income ",
                    ylab = "Higher CDDs ",
                    size = 8)


library(cowplot)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.1, .15, 0.2, 0.2)

ggsave("biv_map.png", finalPlot, height=5, width = 10, scale=2)


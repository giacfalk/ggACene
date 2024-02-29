#############################

# Load packages
library(data.table)
library(plyr)
library(sf)
library(dplyr)
library(FSA)
library(haven)
library(stringr)
library(tidyverse)
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
library(pROC)
library(pbapply)

setwd(wd)

# Load country data
data <- read_rds("supporting_data/global.rds")

###

data <- unique(data[c("country", "adm1")])

###

gadm <- read_sf("supporting_data/gadm_410-levels.gpkg", layer="ADM_1")

gadm <- dplyr::select(gadm, COUNTRY, NAME_1) %>% dplyr::rename(country = COUNTRY, adm1 = NAME_1)

gadm_bk <- gadm
gadm$geom <- NULL

library(fuzzyjoin)

data_c_sp <- stringdist_join(data, gadm, 
                             by = c("country", "adm1"),
                             mode = "left",
                             ignore_case = TRUE, 
                             method = "jw", 
                             max_dist = 99, 
                             distance_col = "dist") %>%
  group_by(country.x, adm1.x) %>%
  slice_min(order_by = data.frame(country.dist, adm1.dist), n = 1)

data_c_sp <- merge(data_c_sp, gadm_bk, by.x=c("country.y", "adm1.y"), by.y=c("country", "adm1"))

gadm <- read_sf("C:/Users/falchetta/OneDrive - IIASA/IIASA_official_RE4AFAGRI_platform/online_dashboards/supporting_files/gadm_410-levels.gpkg", layer="ADM_0")

srs <- ggplot()+
  theme_void()+
  geom_sf(data=st_as_sf(gadm %>% filter(COUNTRY!="Antarctica")), fill="transparent")+
  geom_sf(data=st_as_sf(data_c_sp), fill="orange", colour="black", lwd=0.1)

ggsave("results/graphs_tables/surveyed_regions.png", height = 5, width = 10, scale=1.1, bg="white")

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

## 1) Load libraries and data ##
library(sandwich)
library(lmtest)
library(foreign)
library(ResourceSelection)
library(optmatch)
library(tidyverse)
library(haven)
library(psych)
library(raster)
library(rnaturalearthdata)
library(sf)
library(gdata)
library(exactextractr)
library(nngeo)
library(caret)
library(MatchIt)
library(ggsci)
library(gdata)
library(jtools)
library(glm2)
library(reshape2)
library(cobalt)
library(relaimpo)
library(domir)
library(pscl)
library(margins)

# Set users
user <- 'gf'

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

####

setwd(stub)

load("rscripts/global_spline/results/global_wgt_dmcf.Rdata")


d1 <- reg_ac$data
d2 <- reg_ely$model

d1 <- d1 %>% dplyr::select(country, ac, weight)
d2 <- d2 %>% dplyr::select(country, ln_ely_q, `(weights)`)

d1 <- d1 %>% group_by(country) %>% summarise(ac=weighted.mean(as.numeric(as.character(ac)), weight, na.rm=T))
d2 <- d2 %>% group_by(country) %>% summarise(ely=weighted.mean(exp(ln_ely_q), `(weights)`, na.rm=T))


# 

r <- read_rds("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Kenya/IHBS/kenya_ihbs.rds")
r1 <- r %>% group_by(country) %>% summarise(ac=weighted.mean(as.numeric(as.character(ac)), weight, na.rm=T), ely=weighted.mean(exp(ln_ely_q), weight, na.rm=T))

r1$country <- as.character(r1$country)

dd <- merge(d1, d2, "country")
dd$country <- as.character(dd$country)

dd[25,] <- c("Kenya", r1$ac, r2$ely)

dd <- arrange(dd, as.character(country))

library(xtable)

dd$Surveyyear <- c(2018, 2011, 2018, 2014, 2011, 2014, 2011, 2019, 2017, 2019, 2019, 2019, 2011, 2016, 2020, 2018, 2011, 2014, 2019, 2019, 2011, 2011, 2011, 2018, 2021)

colnames(dd) <- c("Country", "Avg. AC ownership (%)", "Avg ELY cons. (kWh/HH/yr.)", "Survey year")

dd$`Avg. AC ownership (%)` <- round(as.numeric(dd$`Avg. AC ownership (%)`)*100, 2)
dd$`Avg ELY cons. (kWh/HH/yr.)` <- round(as.numeric(dd$`Avg ELY cons. (kWh/HH/yr.)`), 2)

print(xtable(dd, type = "latex"), file = "rscripts/global_spline/results/ac_ely_descriptive_table.tex")




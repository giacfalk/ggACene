rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

## This R-script:
##      1) makes global gridded projection of future AC uptake and utilisation based on the global spline models and on the gridded dataset of expenditure, CDDs/HDDs, and age, gender and sex

## 1) Load libraries and data ##
library(fasterize)
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
library(maptools)
library(countrycode)
library(stars)
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
library(patchwork)

####

setwd(wd)

##############

# extensive margin

###

load(paste0(wd, "/results/xgboost_models_jan24.Rdata"))

###

png("results/graphs_tables/ac_acc_evol.png", height = 600, width = 1200, res = 150)
plot(ac_model)
dev.off()

png("results/graphs_tables/ely_acc_evol.png", height = 600*2, width = 1200*2, res = 300)
plot(ely_model)
dev.off()

###

png("results/graphs_tables/ac_varimp.png", height = 900, width = 1000, res = 220)
plot(varImp(ac_model))
dev.off()

png("results/graphs_tables/ely_varimp.png", height = 900, width = 1000, res = 220)
plot(varImp(ely_model))
dev.off()

rm(list=ls(all=TRUE)) # Removes all previously created variables
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

# Set users
user <- 'fp'
user <- 'gf'
user <- 'gf_server'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

if (user=='gf_server') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

setwd(paste0(stub, "rscripts/global_spline"))

load(paste0(stub, "rscripts/global_spline/results/xgboost_models.Rdata"))

##############

pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions[,2]
}

x <- dplyr::select(global_ac_train_s1, -ac)

# Compute fast (approximate) Shapley values using 10 Monte Carlo repetitions
set.seed(5038)
shap <- fastshap::explain(rrfFit_ac_shapely, X = as.data.frame(x), pred_wrapper = pfun, nsim = 10, shap_only = F, newdata = as.data.frame(x %>% group_by(macroregion) %>% sample_frac(0.01) %>% ungroup())) # test with subsample of newdata
shap_bk <- shap

#########

library(shapviz)

baseline <- mean(pfun(rrfFit_ac_shapely, newdata = global_ac_train_s1)) ## adapt with subsample of newdata
shp <- shapviz(shap_bk, interactions=F, baseline=baseline)

save.image("results/shapely_values_stage1.Rdata")

##########

a <- sv_waterfall(shp, 1:length(shp$X$macroregion))
b <-sv_waterfall(shp, shp$X$macroregion == "Africa")
c <-sv_waterfall(shp, shp$X$macroregion == "Asia")
d <-sv_waterfall(shp, shp$X$macroregion == "Americas")
e <-sv_waterfall(shp, shp$X$macroregion == "Europe")
#f <-sv_waterfall(shp, shp$X$macroregion == "Oceana")

library(patchwork)

a + (b + c + d + e)

ggsave("results/graphs_tables/s1_reg.png", width=10, height=7, scale=1.3)

#

sv_importance(shp, kind = "beeswarm", show_numbers = TRUE)

ggsave("results/graphs_tables/s1.png", width=7, height=7, scale=1)

##############
##############

pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}

x <- dplyr::select(global_ac_train_s2, -ln_ely_q, -ln_ely_q_predicted)

# Compute fast (approximate) Shapley values using 10 Monte Carlo repetitions
set.seed(5038)
shap <- fastshap::explain(rrfFit_ely_shapely, X = as.data.frame(x), pred_wrapper = pfun, nsim = 10, shap_only = F, newdata = as.data.frame(x %>% group_by(macroregion) %>% sample_frac(0.01) %>% ungroup()) ) # test with subsample of newdata
shap_bk <- shap

#########

library(shapviz)

baseline <- mean(pfun(rrfFit_ely_shapely, newdata = global_ac_train_s2)) 
shp <- shapviz(shap_bk, interactions=F, baseline=baseline)

save.image("results/shapely_values_stage2.Rdata")

##########

a <- sv_waterfall(shp, 1:length(shp$X$macroregion))
b <-sv_waterfall(shp, shp$X$macroregion == "Africa")
c <-sv_waterfall(shp, shp$X$macroregion == "Asia")
d <-sv_waterfall(shp, shp$X$macroregion == "Americas")
e <-sv_waterfall(shp, shp$X$macroregion == "Europe")
#f <-sv_waterfall(shp, shp$X$macroregion == "Oceana")

library(patchwork)

a + (b + c + d + e)

ggsave("results/graphs_tables/s2_reg.png", width=10, height=7, scale=1.3)

###

sv_dependence(shp, color_var = "mean_CDD18_db", v="mean_HDD18_db")+
  geom_smooth(method = "lm", formula = y ~ x + I(x^3), size = 1)

#

sv_importance(shp, kind = "beeswarm", show_numbers = TRUE)

ggsave("results/graphs_tables/s2.png", width=7, height=7, scale=1)

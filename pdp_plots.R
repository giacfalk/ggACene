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
library(parallel)
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

load(paste0(stub, "rscripts/global_spline/results/xgboost_models_jan24.Rdata"))

##########################

library(doParallel)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

library(pdp)
 
pdp1 <- pdp::partial(ac_model, 
        pred.var = c("mean_CDD18_db", "ln_total_exp_usd_2011"), 
        trim.outliers = TRUE, chull = FALSE, parallel = TRUE,
        grid.resolution = 30,  prob=T, paropts = list(.packages = "ranger"))

stopCluster(cluster)
gc()

pdp1$mean_CDD18_db <- exp(pdp1$mean_CDD18_db)
pdp1$ln_total_exp_usd_2011 <- exp(pdp1$ln_total_exp_usd_2011)
pdp1$yhat <- (1- pdp1$yhat) *100

png("results/graphs_tables/pdp_1.png", height = 1200, width = 1200, res=200)
pdp::plotPartial(pdp1)
dev.off()

####

pdp2 <- pdp::partial(ely_model, 
             pred.var = c("ln_total_exp_usd_2011", "mean_CDD18_db"), 
             trim.outliers = TRUE, chull = F, parallel = F,
             grid.resolution = 5,  paropts = list(.packages = "ranger"))

pdp2$mean_CDD18_db <- exp(pdp2$mean_CDD18_db)
pdp2$ln_total_exp_usd_2011 <- exp(pdp2$ln_total_exp_usd_2011)
pdp2$yhat <- exp(pdp2$yhat)

pdp3 <- plotPartial(pdp2, levelplot = FALSE, zlab = "AC %", drape = TRUE,
                    colorkey = TRUE, screen = list(z = 20, x = -60))

pdp3

png("results/graphs_tables/pdp_2.png")
pdp3
dev.off()

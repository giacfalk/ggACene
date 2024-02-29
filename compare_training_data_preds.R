# compare training data and prediction data to ensure everything is okay
# i.e., check that projection data are in the training data range

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

## This R-script:
##      1) makes global gridded projection of future AC uptake and utilisation based on the global spline models and on the gridded dataset of expenditure, CDDs/HDDs, and age, gender and sex

## 1) Load libraries and data ##
library(tidyverse)
library(sf)
library(caret)

####

setwd(wd)

##############

# extensive margin

###

load(paste0(wd, "/results/xgboost_models_benchmarks.Rdata"))
load(paste0(wd, "/results/xgboost_models.Rdata"))
load(paste0(wd, "/supporting_data/data_for_global_spline_v2.Rds"))

##########

ssp="SSP2"
year = 2010
rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
orig_data_bk <- shape
output2 <- list()

#

shape$geometry <- NULL
shape <- na.omit(shape)

tapply(exp(ac_model$trainingData$ln_total_exp_usd_2011), ac_model$trainingData$macroregion, summary)
tapply(pull(exp(shape[,paste0("GDP_", ssp, "_", year)])), shape$region, summary)

###################

tapply(ac_model$trainingData$mean_CDD_db, ac_model$trainingData$macroregion, summary)
tapply(log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, 585), "_", ifelse(year==2010, 2015, year))] +1), shape$macroregion, summary)

tapply(ac_model$trainingData$mean_HDD_db, ac_model$trainingData$macroregion, summary)
tapply(log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, 585), "_", ifelse(year==2010, 2015, year))] +1), shape$macroregion, summary)

summary(ac_model$trainingData$edu_head_2)
summary(shape[,paste0("edu_", ssp, "_", year)] )

summary(ac_model$trainingData$age_head)
summary(shape[,paste0("age_", ssp, "_", year)] )

summary(ac_model$trainingData$urban)
summary(shape[,paste0("URB_", ssp, "_", (year))])

#########################

summary(exp(ac_model$trainingData$ln_total_exp_usd_2011))
summary(pull(exp(shape[,paste0("GDP_", ssp, "_", year)])))

a <- ggplot()+
  theme_classic()+
  geom_density(data=ac_model$trainingData, aes(x=ln_total_exp_usd_2011, colour="Training data"), lwd=2)+
  geom_density(data=shape, aes(x=GDP_SSP2_2020, colour="SSP2(45), 2020"))+
  geom_density(data=shape, aes(x=GDP_SSP2_2050, colour="SSP2(45), 2050"))+
  geom_density(data=shape, aes(x=GDP_SSP5_2050, colour="SSP5(85), 2050"))+
  scale_colour_manual(name="", values=c("red", "orange", "blue", "forestgreen"), breaks=c("Training data", "SSP2(45), 2020", "SSP2(45), 2050", "SSP5(85), 2050"))

b <- ggplot()+
  theme_classic()+
  geom_density(data=ac_model$trainingData, aes(x=mean_CDD18_db, colour="Training data"), lwd=2)+
  geom_density(data=shape, aes(x=log(`CDD_245_2020_CMCC-CM2-SR5` + 1), colour="SSP2(45), 2020"))+
  geom_density(data=shape, aes(x=log(`CDD_245_2050_CMCC-CM2-SR5` + 1), colour="SSP2(45), 2050"))+
  geom_density(data=shape, aes(x=log(`CDD_585_2050_CMCC-CM2-SR5` + 1), colour="SSP5(85), 2050"))+
  scale_colour_manual(name="", values=c("red", "orange", "blue", "forestgreen"), breaks=c("Training data", "SSP2(45), 2020", "SSP2(45), 2050", "SSP5(85), 2050"))

c <- ggplot()+
  theme_classic()+
  geom_density(data=ac_model$trainingData, aes(x=mean_HDD18_db, colour="Training data"), lwd=2)+
  geom_density(data=shape, aes(x=log(`HDD_245_2020_CMCC-CM2-SR5` + 1), colour="SSP2(45), 2020"))+
  geom_density(data=shape , aes(x=log(`HDD_245_2050_CMCC-CM2-SR5` + 1), colour="SSP2(45), 2050"))+
  geom_density(data=shape, aes(x=log(`HDD_585_2050_CMCC-CM2-SR5`+1), colour="SSP5(85), 2050"))+
  scale_colour_manual(name="", values=c("red", "orange", "blue", "forestgreen"), breaks=c("Training data", "SSP2(45), 2020", "SSP2(45), 2050", "SSP5(85), 2050"))

d <- ggplot()+
  theme_classic()+
  geom_density(data=ac_model$trainingData, aes(x=urban_sh, colour="Training data"), lwd=2)+
  geom_density(data=shape, aes(x=URB_SSP2_2020, colour="SSP2(45), 2020"))+
  geom_density(data=shape, aes(x=URB_SSP2_2050, colour="SSP2(45), 2050"))+
  geom_density(data=shape, aes(x=URB_SSP5_2050, colour="SSP5(85), 2050"))+
  scale_x_log10()+
  scale_colour_manual(name="", values=c("red", "orange", "blue", "forestgreen"), breaks=c("Training data", "SSP2(45), 2020", "SSP2(45), 2050", "SSP5(85), 2050"))

library(patchwork)

a + b + c + d + plot_annotation(tag_levels = "A", caption="Note: training data include the 25-country pool; SSP data are global and therefore represent all countries") + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal") 

ggsave("results/graphs_tables/ranges_comparison.png", width = 10, height = 8, scale=1)

###

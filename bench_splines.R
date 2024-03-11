#############################
# random forests

training_frac = 1

# Load packages
library(data.table)
library(plyr)
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

#

setwd(wd)

# Load country data
load("results/global_wgt_dmcf.RData")

global = reg_ely$data
global <- as.data.frame(global)

###

l <- list.files(path="hurs", pattern="rds", full.names = T)
l <- lapply(l, read_rds)

for(i in 1:length(l)){
  print(i)
  l[[i]] <- dplyr::select(l[[i]], state, country, hurs)
}

l <- dplyr::bind_rows(l)
l <- group_by(l, country, state) %>% dplyr::summarise(hurs=mean(hurs,na.rm=T))

global_bk <- global
global$country <- as.character(global$country)
global$id <- 1:nrow(global)
global <- merge(global, l, by.x=c("country", "adm1"), by.y=c("country", "state"))
global <- global[!duplicated(global$id),]
global$id <- NULL

###

not_all_na <- function(x) any(!is.na(x))
not_any_na <- function(x) all(!is.na(x))

global <- dplyr::select(global, country, ac, ln_ely_q, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, n_members, ely_p_usd_2011, weight, hurs)

global <- dplyr::select(global, where(not_all_na))

gc()

######################

# adjustment expenditure to gdp per capita

adj_gtap <- read.csv("adj_gtap/exp_gdp_capita_markups.csv")

adj_gtap$REG <- toupper(adj_gtap$REG)

global$iso3c <- countrycode::countrycode(global$country, 'country.name', 'iso3c')

global <- merge(global, adj_gtap, by.x="iso3c", by.y="REG")

#global$ln_total_exp_usd_2011 <- log(exp(global$ln_total_exp_usd_2011) / global$markup)

# add a macroregion variable

setwd(wd)

load(paste0(wd, "/supporting_data/adj_factors.Rds"))

ss <- dplyr::select(ss, country, region)
colnames(ss)[2] <- "macroregion"

global <- merge(global, ss, by="country")

global$macroregion <- countrycode::countrycode(global$country, 'country.name', 'continent')

#

# remove outliers

# global <- filter(global, ln_ely_q<10 & ln_ely_q > 0)
# global <- filter(global, ln_total_exp_usd_2011<12.5& ln_total_exp_usd_2011 > 0)
# global <- filter(global, n_members<10)
# global <- filter(global, age_head<99)

global <- global %>% group_by(country) %>% filter_if(is.numeric, all_vars(between(., quantile(., .01, na.rm=T), quantile(., .99, na.rm=T))))

global$mean_CDD18_db <- log((global$mean_CDD18_db*100)+1)
global$mean_HDD18_db <- log((global$mean_HDD18_db*100)+1)

# convert some factors to numeric to allow projections
global$edu_head_2 <- as.numeric(as.character(global$edu_head_2))

##############
library(sf)
ctrs <- st_as_sf(rnaturalearthdata::countries110)
ctrs <- st_transform(ctrs, 3395)
ctrs <- st_make_valid(ctrs)
ctrs$area <- as.numeric(st_area(ctrs)/1e6)
ctrs <- st_transform(ctrs, 4326)
ctrs <- st_make_valid(ctrs)
ctrs$X <- as.vector(st_coordinates(st_centroid(ctrs))[,1])
ctrs$Y <- as.vector(st_coordinates(st_centroid(ctrs))[,2])
ctrs <- dplyr::select(ctrs, iso_a3, area, X, Y)
ctrs$geometry <- NULL

global$X <- NULL
global <- merge(global, ctrs, by.x="iso3c", by.y="iso_a3")

global_ac <- global

####

terc_inc <- quantile(exp(global_ac$ln_total_exp_usd_2011), seq(0, 1, 0.33))
terc_cdd <- quantile(exp(global_ac$mean_CDD18_db), seq(0, 1, 0.33))
terc_hdd <- quantile(exp(global_ac$mean_HDD18_db), seq(0, 1, 0.33))

global_ac$inc_q <- cut(exp(global_ac$ln_total_exp_usd_2011), terc_inc, include.lowest = T, labels=c(1,2,3))
global_ac$cdd_q <- cut(exp(global_ac$mean_CDD18_db), terc_cdd, include.lowest = T, labels=c(1,2,3))
global_ac$hdd_q <-  cut(exp(global_ac$mean_HDD18_db), terc_hdd, include.lowest = T, labels=c(1,2,3))

global_ac$hdd_q <- ifelse(exp(global_ac$mean_HDD18_db)<terc_hdd[1], 1, global_ac$hdd_q)
global_ac$hdd_q <- ifelse(exp(global_ac$mean_HDD18_db)>terc_hdd[4], 3, global_ac$hdd_q)
global_ac$hdd_q <- as.factor(global_ac$hdd_q)

global_ac$cdd_q <- ifelse(exp(global_ac$mean_CDD18_db)<terc_cdd[1], 1, global_ac$cdd_q)
global_ac$cdd_q <- ifelse(exp(global_ac$mean_CDD18_db)>terc_cdd[4], 3, global_ac$cdd_q)
global_ac$cdd_q <- as.factor(global_ac$cdd_q)

global_ac$inc_q <- ifelse(exp(global_ac$ln_total_exp_usd_2011)<terc_inc[1], 1, global_ac$inc_q)
global_ac$inc_q <- ifelse(exp(global_ac$ln_total_exp_usd_2011)>terc_inc[4], 3, global_ac$inc_q)
global_ac$inc_q <- as.factor(global_ac$inc_q)

#######

#global_ac %>% group_by(country) %>% summary()

######

# sample (for demo quick runs)
global_ac <- global_ac %>% group_by(macroregion, inc_q, cdd_q, hdd_q) %>% sample_frac(training_frac) %>% ungroup()

#global_ac <- filter(global_ac, inc_q==3)

#global_ac$ln_total_exp_usd_2011 = exp(global_ac$ln_total_exp_usd_2011)

# country-stratified random sampling

library(splitstackshape)
set.seed(2022)

global_ac_train <- stratified(global_ac, c("country", "macroregion", "inc_q", "cdd_q", "hdd_q"), 0.7, bothSets = T)

#global_ac_train <- stratified(global_ac, c("ac", "inc_q", "cdd_q", "hdd_q"), 0.7, bothSets = T)

global_ac_test <- as.data.frame(global_ac_train[[2]])
global_ac_train <- as.data.frame(global_ac_train[[1]])

#semiArtificial::dataSimilarity(global_ac_train %>% dplyr::select(where(is.numeric)), global_ac_test %>% dplyr::select(where(is.numeric)))

global_ac_train$ac<-gsub('0','No',global_ac_train$ac)
global_ac_train$ac<-gsub('1','Yes',global_ac_train$ac)

global_ac_test$ac<-gsub('0','No',global_ac_test$ac)
global_ac_test$ac<-gsub('1','Yes',global_ac_test$ac)

global_ac_train$ac <- as.factor(global_ac_train$ac)
global_ac_test$ac <- as.factor(global_ac_test$ac)

levels(global_ac_train$ac) <- make.names(levels(factor(global_ac_train$ac)))
levels(global_ac_test$ac) <- make.names(levels(factor(global_ac_test$ac)))

prop.table(table(global_ac_train$ac))
prop.table(table(global_ac_test$ac))

#

library(caret)
library(parallel)
library(doParallel)

global_ac_train_s1 <- dplyr::select(global_ac_train, ac, macroregion, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, ely_p_usd_2011, weight, hurs)

#global_ac_train_s1 <- dplyr::select(global_ac_train, ac, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, inc_q, cdd_q, hdd_q, weight)

global_ac_test_s1 <- dplyr::select(global_ac_test, ac, macroregion, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2, age_head, ely_p_usd_2011, weight, hurs)


#model_weights <- global_ac_train$weight


##################################

#        Extensive margin        #

##################################

# Enable multithread support
# cores_2_use <- floor(0.65*detectCores())
# cl <- makeCluster(cores_2_use, outfile = "parallel_log4.txt")
# registerDoParallel(cl)

###

#global_ac_train_s1$strata <- paste(global_ac_train_s1$inc_q, global_ac_train_s1$cdd_q, global_ac_train_s1$hdd_q, sep="_")

library(CAST)

# folds <- 5
# cvIndex <- CAST::CreateSpacetimeFolds(global_ac_train_s1, spacevar="strata", k=folds)

# country_idx <- global_ac_train_s1$country
# weight_idx <- global_ac_train_s1$weight

global_ac_train_s1$strata <- NULL
global_ac_train_s1$country <- NULL
#global_ac_train_s1$weight <- NULL

library(mltools)

# MCC <- function(data, lev = NULL, model = NULL) {
#   
#   pred_ac <- data.frame(global_ac_train_s1)
#   
#   pred_ac$ac <- ifelse(data$pred==lev[1], 1,0)
#   
#   pred_ac$country <- country_idx
#   pred_ac$weight <- weight_idx
#   
#   obs_ac <- data.frame(global_ac_train_s1)
#   
#   obs_ac$ac <- ifelse(data$obs==lev[1], 1,0)
#   
#   obs_ac$country <- country_idx
#   obs_ac$weight <- weight_idx
#   
#   pred_ac <- pred_ac %>% group_by(country) %>% dplyr::summarise(ac=weighted.mean(ac, weight, na.rm=T)) %>% ungroup() %>% pull(ac)
#   
#   obs_ac <- obs_ac %>% group_by(country) %>% dplyr::summarise(ac=weighted.mean(ac, weight, na.rm=T)) %>% ungroup() %>% pull(ac)
#   
#   mse_ac <- mean((pred_ac-obs_ac)^2, na.rm=T)
#   
#   c(MCC = mse_ac)
# }

folds <- 10

# Model cross validations options
fitControl <- trainControl(
  ## Repeated 5-fold CV 
  method = "cv",
  number = folds,
  #index = cvIndex$index,
  #indexOut = cvIndex$indexOut,
  verboseIter = TRUE, allowParallel = F, classProbs = T)


model_weights <- ifelse(global_ac_train_s1$ac == "1",
                        (1/table(global_ac_train_s1$ac)[1]) * 0.5,
                        (1/table(global_ac_train_s1$ac)[2]) * 0.5)

#fitControl$sampling <- "smote"

# https://stats.stackexchange.com/questions/301364/gam-optimization-methods-in-mgcv-r-package-which-to-choose

# definisci manualmente i parameter (num.trees, mtry, splitrule) space con tunematrix
rrfFit_ac_glm <- caret::train(ac ~ .,
                data = as.data.frame(global_ac_train_s1),
                method = 'glm',
                weights = model_weights,
                metric="Kappa",
                trControl = fitControl,
                na.action = na.omit,
                maximize = T)

print(rrfFit_ac_glm)

##############

tgrid <- expand.grid(
  select = c("TRUE", "FALSE"),
  method = c("GCV.Cp")
)

# definisci manualmente i parameter (num.trees, mtry, splitrule) space con tunematrix
rrfFit_ac_gam <- caret::train(ac ~ .,
                data = as.data.frame(global_ac_train_s1),
                method = 'gam',
                importance = "permutation",
                tuneGrid = tgrid,
                weights = model_weights,
                metric="Kappa",
                trControl = fitControl,
                na.action = na.omit,
                maximize = T)



# Print CV accuracy
print(rrfFit_ac_gam)

#######

library(rpart)
library(rpart.plot)

fit.tree = rpart(ac ~ ., data=as.data.frame(global_ac_train_s1), method = "class", cp=0.008)

png("graphs_tables/cart_tree.png", height = 1200, width = 1000)
rpart.plot(fit.tree)
dev.off()

#################

ac_model_glm <- rrfFit_ac_glm
ac_model_gam <- rrfFit_ac_gam

#########################

#as.numeric(prCalibrate(ifelse(global_ac_train_s1$ac=="Yes", 1, 0), predict(ac_model, newdata=global_ac_train_s1, type="prob")[,2])$cal.probs)

# Predicted probabilities
global_ac_train$phat0_obs <- predict(ac_model_glm, newdata=global_ac_train_s1, type="prob")[,2]

global_ac_test$phat0_obs <- predict(ac_model_glm, newdata=global_ac_test_s1, type="prob")[,2]

####################################################

global_ac_train_s2 <- dplyr::select(global_ac_train, country, macroregion, ln_ely_q, phat0_obs, mean_CDD18_db,  mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2, age_head, ely_p_usd_2011, weight, hurs)

folds <- 10


# Model cross validations options
fitControl <- trainControl(
  ## Repeated 5-fold CV 
  method = "cv",
  number = folds,
  # index = cvIndex$index,
  # indexOut = cvIndex$indexOut,
  verboseIter = TRUE,
  returnResamp = "all", allowParallel = F)

global_ac_train_s2$country <- NULL

model_weights <- global_ac_train$weight

# definisci manualmente i parameter (num.trees, mtry, splitrule) space con tunematrix
rrfFit_ely_lm <- caret::train(ln_ely_q ~ .,
                    data = as.data.frame(global_ac_train_s2),
                    method = 'lm',
                    metric = "Rsquared",
                    trControl = fitControl,
                    weights = model_weights,
                    importance = "permutation",
                    na.action = na.omit)


################

global_ac_train$phat0_obs <- predict(ac_model_gam, newdata=global_ac_train_s1, type="prob")[,2]
global_ac_test$phat0_obs <- predict(ac_model_gam, newdata=global_ac_test_s1, type="prob")[,2]

################

global_ac_train_s2 <- dplyr::select(global_ac_train, country, macroregion, ln_ely_q, phat0_obs, mean_CDD18_db,  mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2, age_head, ely_p_usd_2011, weight, hurs)

tgrid <- expand.grid(
  select = c("TRUE", "FALSE"),
  method = c("GCV.Cp")
)

global_ac_train_s2$country <- NULL

rrfFit_ely_gam <- caret::train(ln_ely_q ~ .,
                    data = as.data.frame(global_ac_train_s2),
                    method = 'gam',
                    metric = "Rsquared",
                    tuneGrid = tgrid,
                    trControl = fitControl,
                    weights = model_weights,
                    importance = "permutation",
                    na.action = na.omit)


##########################
##########################

training_acc_ac_glm <- auc(roc(factor(global_ac_train_s1$ac, ordered = T), factor(predict(rrfFit_ac_glm, global_ac_train_s1), ordered = T)))
testing_acc_ac_glm <- auc(roc(factor(global_ac_test_s1$ac, ordered = T), factor(predict(rrfFit_ac_glm, global_ac_test_s1), ordered = T)))

training_kappa_ac_glm <- psych::cohen.kappa(table(factor(predict(rrfFit_ac_glm, global_ac_train_s1), ordered = T), factor(global_ac_train_s1$ac, ordered = T)))$kappa
testing_kappa_ac_glm <- psych::cohen.kappa(table(factor(predict(rrfFit_ac_glm, global_ac_test_s1), ordered = T), factor(global_ac_test_s1$ac, ordered = T)))$kappa

training_acc_ac_gam <- auc(roc(factor(global_ac_train_s1$ac, ordered = T), factor(predict(rrfFit_ac_gam, global_ac_train_s1), ordered = T)))
testing_acc_ac_gam <- auc(roc(factor(global_ac_test_s1$ac, ordered = T), factor(predict(rrfFit_ac_gam, global_ac_test_s1), ordered = T)))

training_kappa_ac_gam <- psych::cohen.kappa(table(factor(predict(rrfFit_ac_gam, global_ac_train_s1), ordered = T), factor(global_ac_train_s1$ac, ordered = T)))$kappa
testing_kappa_ac_gam <- psych::cohen.kappa(table(factor(predict(rrfFit_ac_gam, global_ac_test_s1), ordered = T), factor(global_ac_test_s1$ac, ordered = T)))$kappa

setwd(wd)
load("results/xgboost_models_jan24.Rdata")
load("results/xgboost_models_benchmarks_jan24.Rdata")

training_acc_ac_rf <- auc(roc(factor(ac_model$trainingData$.outcome, ordered = T), factor(predict(ac_model, ac_model$trainingData), ordered = T)))
testing_acc_ac_rf <- auc(roc(factor(global_ac_test_s1$ac, ordered = T), factor(predict(ac_model, global_ac_test_s1), ordered = T)))

training_kappa_ac_rf <- psych::cohen.kappa(table(factor(predict(ac_model, ac_model$trainingData), ordered = T), factor(ac_model$trainingData$.outcome, ordered = T)))$kappa
testing_kappa_ac_rf <- max(ac_model$results$Kappa)
  
##############

library(modelsummary)

datafr <-  expand.grid(Model=c("GLM", "GAM", "RF"), Set=c("Training", "Testing"), metric=c("Kappa", "AUC"))
datafr$value <- c(training_kappa_ac_glm, training_kappa_ac_gam, training_kappa_ac_rf, testing_kappa_ac_glm, testing_kappa_ac_gam, testing_kappa_ac_rf, training_acc_ac_glm, training_acc_ac_gam, training_acc_ac_rf, testing_acc_ac_glm, testing_acc_ac_gam, testing_acc_ac_rf)

datasummary(Model * Set ~ metric * value * mean, data = datafr, output = "results/graphs_tables/bench_table_1.tex")

#############

mse <- function(real,modeled){
  round(mean((real - modeled)^2, na.rm=T), 2)
}


training_r2_ely_glm <- cor(global_ac_train_s2$ln_ely_q, predict(rrfFit_ely_lm, global_ac_train_s2))^2
testing_r2_ely_glm <- cor(global_ac_test$ln_ely_q, predict(rrfFit_ely_lm, global_ac_test))^2

training_mse_ely_glm <- mse(global_ac_train_s2$ln_ely_q, predict(rrfFit_ely_lm, global_ac_train_s2))
testing_mse_ely_glm <- mse(global_ac_test$ln_ely_q, predict(rrfFit_ely_lm, global_ac_test))

training_r2_ely_gam <- cor(global_ac_train_s2$ln_ely_q, predict(rrfFit_ely_gam, global_ac_train_s2))^2
testing_r2_ely_gam <- cor(global_ac_test$ln_ely_q, predict(rrfFit_ely_gam, global_ac_test))^2

training_mse_ely_gam <- mse(global_ac_train_s2$ln_ely_q, predict(rrfFit_ely_gam, global_ac_train_s2))
testing_mse_ely_gam <- mse(global_ac_test$ln_ely_q, predict(rrfFit_ely_gam, global_ac_test))

training_r2_ely_rf <- cor(ely_model$trainingData$.outcome, predict(ely_model, ely_model$trainingData))^2
testing_r2_ely_rf <- max(ely_model$results$Rsquared)

training_mse_ely_rf <-  mse(ely_model$trainingData$.outcome, predict(ely_model, ely_model$trainingData))
testing_mse_ely_rf <- ely_model$finalModel$prediction.error

####
  
datafr <-  expand.grid(Model=c("LM", "GAM", "RF"), Set=c("Training", "Testing"), metric=c("R-squared", "MSE"))
datafr$value <- c(training_r2_ely_glm, training_r2_ely_gam, training_r2_ely_rf, testing_r2_ely_glm, testing_r2_ely_gam, testing_r2_ely_rf, training_mse_ely_glm, training_mse_ely_gam, training_mse_ely_rf, testing_mse_ely_glm, testing_mse_ely_gam, testing_mse_ely_rf)

datasummary(Model * Set ~ metric * value * mean, data = datafr, output = "results/graphs_tables/bench_table_2.tex")

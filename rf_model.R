#############################
# random forests

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

tapply(setdiff(global_bk, global %>% dplyr::select(-hurs))$adm1, setdiff(global_bk, global %>% dplyr::select(-hurs))$country, unique)

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

setwd(wd)

load(paste0(wd, "supporting_data/adj_factors.Rds"))

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

global <- na.omit(global)

global <- global %>% group_by(country) %>% filter_if(is.numeric, all_vars(between(., quantile(., .01), quantile(., .99))))

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

global_ac_train_s1 <- dplyr::select(global_ac_train, ac, macroregion, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, weight, ely_p_usd_2011, hurs)

#global_ac_train_s1 <- dplyr::select(global_ac_train, ac, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, inc_q, cdd_q, hdd_q, weight)

global_ac_test_s1 <- dplyr::select(global_ac_test, ac, macroregion, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2, age_head, weight, ely_p_usd_2011, hurs)


#model_weights <- global_ac_train$weight


##################################

#        Extensive margin        #

##################################

# Enable multithread support
# cores_2_use <- floor(0.65*detectCores())
# cl <- makeCluster(cores_2_use, outfile = "parallel_log3.txt")
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

tgrid <- expand.grid(
  mtry = 2:8,
  splitrule = c("gini", "extratrees"),
  min.node.size = seq(1, 20, 1)
)

# definisci manualmente i parameter (num.trees, mtry, splitrule) space con tunematrix
rrfFit <- caret::train(ac ~ .,
                data = as.data.frame(global_ac_train_s1),
                method = 'ranger',
                num.trees  = 100,
                importance = "permutation",
                tuneGrid = tgrid,
                weights = model_weights,
                metric="Kappa",
                trControl = fitControl,
                na.action = na.omit,
                maximize = T)

library(ranger)

rrfFit_ac_shapely <- ranger(ac ~ .,
                            data = global_ac_train_s1, probability = T)


# rrfFit <- train(ac ~ ., 
#                 data = as.data.frame(global_ac_train_s1),
#                 method = 'xgbTree',
#                 tuneLength = 5,
#                 weights = model_weights,
#                 metric="Kappa",
#                 trControl = fitControl,
#                 na.action = na.omit,
#                 maximize = T)

# Print CV accuracy
print(rrfFit)

psych::cohen.kappa(table(factor(predict(rrfFit, global_ac_train_s1), ordered = T), factor(global_ac_train_s1$ac, ordered = T)))$kappa

# # country level accuracy
# for (mr in unique(global_ac_train_s1$macroregion)){
#   
#   print(paste(mr, psych::cohen.kappa(table(factor(predict(rrfFit, global_ac_train_s1[global_ac_train_s1$macroregion==mr,]), ordered = T), factor(global_ac_train_s1$ac[global_ac_train_s1$macroregion==mr], ordered = T)))$kappa))
# }


##########

roc_obj <- roc(factor(global_ac_train_s1$ac, ordered = T), factor(predict(rrfFit, newdata=global_ac_train_s1), ordered = T))
plot(roc_obj, col="red", lwd=3, main="ROC curve QDA")
auc_out <-   auc(roc_obj)
auc_out
training_acc_ac <- auc_out

# Variable importance
varImp_tr <- varImp(rrfFit)
plot(varImp_tr)

# Predict on test
global_ac_test_s1$ac_predicted = (predict(rrfFit, global_ac_test_s1, type="prob"))[,2]

outer <- list()

for (thresh in seq(0.01, 0.99, 0.01)){

outer[[as.character(thresh)]] <- psych::cohen.kappa(table(factor(ifelse(global_ac_test_s1$ac_predicted>thresh, "Yes", "No"), ordered = T), factor(global_ac_test_s1$ac, ordered = T)))$kappa

top_threshold <- as.numeric(names(which.max(unlist(outer))))

}
  
# AUC for test (= test accuracy)
roc_obj <- roc(factor(global_ac_test_s1$ac, ordered = T), factor(ifelse(global_ac_test_s1$ac_predicted>top_threshold, "Yes", "No"), ordered = T))
plot(roc_obj, col="red", lwd=3, main="ROC curve QDA")
auc_out <-   auc(roc_obj)
auc_out
testing_acc_ac <- auc_out

psych::cohen.kappa(table(factor(ifelse(global_ac_test_s1$ac_predicted>top_threshold, "Yes", "No"), ordered = T), factor(global_ac_test_s1$ac, ordered = T)))$kappa

###########


# for (mr in unique(global_ac_train_s1$macroregion)){
#   
#   print(paste(mr, psych::cohen.kappa(table(factor(predict(rrfFit, global_ac_test_s1[global_ac_test_s1$macroregion==mr,]), ordered = T), factor(global_ac_test_s1$ac[global_ac_test_s1$macroregion==mr], ordered = T)))$kappa))
# }

# Shut down parallel
# stopCluster(cl)
# registerDoSEQ()

######################

#select probability threshold that maximises cohen's kappa'

# out <- list()
# 
# for (thresh in seq(0.01, 0.9, 0.01)){
#   
#   out[[as.character(thresh)]] <- psych::cohen.kappa(table(factor(global_ac_test_s1$ac, ordered = T), factor(ifelse((global_ac_test_s1$ac_predicted >thresh)==TRUE, "Yes", "No"), ordered = T)))$kappa
#   
# }

#source("prCalibrate.R")

ac_model <- rrfFit

#########################

#as.numeric(prCalibrate(ifelse(global_ac_train_s1$ac=="Yes", 1, 0), predict(ac_model, newdata=global_ac_train_s1, type="prob")[,2])$cal.probs)

# Predicted probabilities
global_ac_train$phat0_obs <- predict(ac_model, newdata=global_ac_train_s1, type="prob")[,2]

# # Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
# global_ac_train$ac_obs <- ifelse(global_ac_train$phat0_obs>0.5 & !is.na(global_ac_train$phat0_obs), 1 , 0)
# #table(global_ac_train$ac, global_ac_train$ac_obs)
# 
# # Selection term for intensive margin part
# global_ac_train$xb_noac = 1-global_ac_train$phat0_obs               
# global_ac_train$selection = ifelse(global_ac_train$ac==1, 
#                                    (global_ac_train$xb_noac*log(global_ac_train$xb_noac)/global_ac_train$phat0_obs) + log(global_ac_train$phat0_obs), 
#                                    (global_ac_train$phat0_obs*log(global_ac_train$phat0_obs)/global_ac_train$xb_noac) + log(global_ac_train$xb_noac))
# 
# 
# #
# 
global_ac_test$phat0_obs <- predict(ac_model, newdata=global_ac_test_s1, type="prob")[,2]
# 
# 
# # Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
# global_ac_test$ac_obs <- ifelse(global_ac_test$phat0_obs>0.5 & !is.na(global_ac_test$phat0_obs), 1 , 0)
# #table(global_ac_test$ac, global_ac_test$ac_obs)
# 
# # Selection term for intensive margin part
# global_ac_test$xb_noac = 1-global_ac_test$phat0_obs               
# global_ac_test$selection = ifelse(global_ac_test$ac==1, 
#                                    (global_ac_test$xb_noac*log(global_ac_test$xb_noac)/global_ac_test$phat0_obs) + log(global_ac_test$phat0_obs), 
#                                    (global_ac_test$phat0_obs*log(global_ac_test$phat0_obs)/global_ac_test$xb_noac) + log(global_ac_test$xb_noac))
# 
rm(rrfFit)

####################################################
####################################################

#####################

# Enable multithread support
# cores_2_use <- floor(0.5*detectCores())
# cl <- makeCluster(cores_2_use, outfile = "parallel_log_2.txt")
# registerDoParallel(cl)

global_ac_train_s2 <- dplyr::select(global_ac_train, country, macroregion, ln_ely_q, phat0_obs, mean_CDD18_db,  mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2, age_head, weight, ely_p_usd_2011, hurs)

# folds <- 20
# cvIndex <- CAST::CreateSpacetimeFolds(global_ac_train_s2, spacevar="country", k=folds, class="ln_ely_q")

folds <- 10

tgrid <- expand.grid(
  mtry = 2:9,
  splitrule="extratrees",
  min.node.size = seq(1, 20, 1)
)

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
rrfFit_ely <- caret::train(ln_ely_q ~ .,
                    data = as.data.frame(global_ac_train_s2),
                    method = 'ranger',
                    metric = "Rsquared",
                    num.trees = 500,
                    tuneGrid = tgrid,
                    trControl = fitControl,
                    weights = model_weights,
                    importance = "permutation",
                    na.action = na.omit)

rrfFit_ely_shapely <- ranger(ln_ely_q ~ .,
                             data = global_ac_train_s2)

# rrfFit_ely <- train(ln_ely_q ~ ., 
#                     data = as.data.frame(global_ac_train_s2),
#                     method = 'xgbTree',
#                     metric = "Rsquared",
#                     tuneLength = 5,
#                     trControl = fitControl,
#                     weights = model_weights,
#                     importance = "permutation",
#                     na.action = na.omit)

# Print CV accuracy
print(rrfFit_ely)

##########

r2_train <- cor(global_ac_train$ln_ely_q, predict(rrfFit_ely, global_ac_train_s2))^2
r2_train

# Variable importance
varImp_test <- varImp(rrfFit_ely)

# Predict on test
global_ac_test$ln_ely_q_predicted = (predict(rrfFit_ely, global_ac_test))
r2_test <- cor(global_ac_test$ln_ely_q, global_ac_test$ln_ely_q_predicted)^2
r2_test

# Shut down parallel
# stopCluster(cl)
# registerDoSEQ()

######################

ely_model <- rrfFit_ely
rm(rrfFit_ely)

training_acc_ac <- as.numeric(training_acc_ac)
testing_acc_ac <- as.numeric(testing_acc_ac)

save("training_acc_ac", "testing_acc_ac", "varImp_tr", "varImp_test", "r2_train", "r2_test", file=paste0(wd, "results/xgboost_models_benchmarks_jan24.Rdata"))

#

global = bind_rows(global_ac_train, global_ac_test)  

#

save("ac_model", "ely_model", "rrfFit_ac_shapely", "rrfFit_ely_shapely", "global_ac_train_s1", "global_ac_train_s2", "global", "global_ac", "terc_inc", "terc_cdd", "terc_hdd", file=paste0(wd, "results/xgboost_models_jan24.Rdata"))

###

mape <- function(real,modeled){
  
  round(mean(abs((real - modeled)/real), na.rm=T), 2)
  
}


mse <- function(real,modeled){
  
  round(mean((real - modeled)^2, na.rm=T), 2)
  
}


global_ac_train_s2$ln_ely_q_predicted <- (predict(ely_model, global_ac_train_s2))

val_2_tr <- ggplot(global_ac_train_s2, aes(x=ln_ely_q_predicted, y=ln_ely_q)) + geom_abline() + geom_point(pch='.') + ggtitle("Training set, 2nd stage (intensive margin)")+
  annotate("text", x=5.5, y=10, label= paste0("r2 = ", round(cor(global_ac_train_s2$ln_ely_q_predicted, global_ac_train_s2$ln_ely_q)^2, 2)))+
  annotate("text", x=5.5, y=9.5, label= paste0("MSE = ", mse(global_ac_train_s2$ln_ely_q_predicted, global_ac_train_s2$ln_ely_q)))+
  xlab("Predicted electricity consumption")+
  ylab("Observed electricity consumption")

val_2_te <- ggplot(global_ac_test, aes(x=ln_ely_q_predicted, y=ln_ely_q)) + geom_abline() + geom_point(pch='.') + ggtitle("Testing set, 2nd stage (intensive margin)")+
  annotate("text", x=5.5, y=10, label= paste0("r2 = ", round(cor(global_ac_test$ln_ely_q_predicted, global_ac_test$ln_ely_q)^2, 2)))+
  annotate("text", x=5.5, y=9.5, label= paste0("MSE = ", mse(global_ac_test$ln_ely_q_predicted, global_ac_test$ln_ely_q)))+
  xlab("Predicted electricity consumption")+
  ylab("Observed electricity consumption")

library(patchwork)
library(scales)

cfm_tr <- confusionMatrix(predict(ac_model, global_ac_train_s1), global_ac_train_s1$ac)
cfm_te <- confusionMatrix(predict(ac_model, global_ac_test_s1), global_ac_test_s1$ac)

cfm_tr$overall[8] <- training_acc_ac
cfm_te$overall[8] <- testing_acc_ac

ggplotConfusionMatrix1 <- function(m){
  mytitle <- paste0("Accuracy ", percent_format()(m$overall[1]), ", ",
                    "Kappa ", percent_format()(m$overall[2]), ", ",
                    "AUC ", percent_format()(m$overall[8]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") +
    ggtitle("Training set, 1st stage (extensive margin)", subtitle = mytitle)
  return(p)
}

ggplotConfusionMatrix2 <- function(m){
  mytitle <- paste0("Accuracy ", percent_format()(m$overall[1]), ", ",
                    "Kappa ", percent_format()(m$overall[2]), ", ",
                    "AUC ", percent_format()(m$overall[8]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") +
    ggtitle("Testing set, 1st stage (extensive margin)", subtitle = mytitle)
  return(p)
}


val_1_tr <- ggplotConfusionMatrix1(cfm_tr)
val_1_te <- ggplotConfusionMatrix2(cfm_te)

###

val_1_tr + val_1_te + val_2_tr + val_2_te + plot_layout(ncol=2)

ggsave("results/graphs_tables/accuracy_plot_jan24.png", scale=1.5, height = 5, width = 7)

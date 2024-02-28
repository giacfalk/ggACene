rm(list=ls(all=TRUE)) # Removes all previously created variables
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

# Set users
user <- 'fp'
user <- 'gf'
user <- 'gf_server'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf_server') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

setwd(paste0(stub, "6-Projections/rscripts/global_spline"))

# descriptive statistics
#######

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

##
# calculate IQR of GCMs

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ac <- bind_cols(shape_ac, future_acc_q1, future_acc_q3)

#####

shape_ely_diff <- shape_ely_diff %>% dplyr::select(id, contains(".20") & contains("SSP"))
shape_ely_diff$region <- NULL
colnames(shape_ely_diff)[2:261] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:261])

shape_ely_diff <- ungroup(shape_ely_diff)

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ely_diff <- bind_cols(shape_ely_diff, future_acc_q1, future_acc_q3)

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff, "id")
shape_ely_diff_m[is.na(shape_ely_diff_m)] = 0

########################

# x axis: percentage change in total energy use (from negative to >100)
# y axis: million people (cumulative sum)
# facets: region
# lines coloured by SSP/decile

shape_ely_diff_m$cons_AC_SSP2.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020), 0, shape_ely_diff_m$cons_AC_SSP2.2020)
shape_ely_diff_m$cons_AC_SSP2.2050 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2050), 0, shape_ely_diff_m$cons_AC_SSP2.2050)
shape_ely_diff_m$pct_chg_cons_AC_SSP2 = (shape_ely_diff_m$cons_AC_SSP2.2050 - shape_ely_diff_m$cons_AC_SSP2.2020) / shape_ely_diff_m$cons_AC_SSP2.2020
shape_ely_diff_m$pct_chg_cons_AC_SSP2 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP2) | shape_ely_diff_m$pct_chg_cons_AC_SSP2 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP2)
shape_ely_diff_m$pct_chg_cons_AC_SSP2 = ifelse(shape_ely_diff_m$cons_AC_SSP2.2050==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP2)
shape_ely_diff_m$pct_chg_cons_AC_SSP2 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP2<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP2)

shape_ely_diff_m$cons_AC_SSP2.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP2.2020_q1)
shape_ely_diff_m$cons_AC_SSP2.2050_q1 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2050_q1), 0, shape_ely_diff_m$cons_AC_SSP2.2050_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 = (shape_ely_diff_m$cons_AC_SSP2.2050_q1 - shape_ely_diff_m$cons_AC_SSP2.2020_q1) / shape_ely_diff_m$cons_AC_SSP2.2020_q1
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1) | shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 = ifelse(shape_ely_diff_m$cons_AC_SSP2.2050_q1==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1)

shape_ely_diff_m$cons_AC_SSP2.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP2.2020_q3)
shape_ely_diff_m$cons_AC_SSP2.2050_q3 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2050_q3), 0, shape_ely_diff_m$cons_AC_SSP2.2050_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 = (shape_ely_diff_m$cons_AC_SSP2.2050_q3 - shape_ely_diff_m$cons_AC_SSP2.2020_q3) / shape_ely_diff_m$cons_AC_SSP2.2020_q3
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3) | shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 = ifelse(shape_ely_diff_m$cons_AC_SSP2.2050_q3==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3)

shape_ely_diff_m$cons_AC_SSP5.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020), 0, shape_ely_diff_m$cons_AC_SSP5.2020)
shape_ely_diff_m$cons_AC_SSP5.2050 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2050), 0, shape_ely_diff_m$cons_AC_SSP5.2050)
shape_ely_diff_m$pct_chg_cons_AC_SSP5 = (shape_ely_diff_m$cons_AC_SSP5.2050 - shape_ely_diff_m$cons_AC_SSP5.2020) / shape_ely_diff_m$cons_AC_SSP5.2020
shape_ely_diff_m$pct_chg_cons_AC_SSP5 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP5) | shape_ely_diff_m$pct_chg_cons_AC_SSP5 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP5)
shape_ely_diff_m$pct_chg_cons_AC_SSP5 = ifelse(shape_ely_diff_m$cons_AC_SSP5.2050==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP5)
shape_ely_diff_m$pct_chg_cons_AC_SSP5 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP5<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP5)

shape_ely_diff_m$cons_AC_SSP5.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP5.2020_q1)
shape_ely_diff_m$cons_AC_SSP5.2050_q1 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2050_q1), 0, shape_ely_diff_m$cons_AC_SSP5.2050_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 = (shape_ely_diff_m$cons_AC_SSP5.2050_q1 - shape_ely_diff_m$cons_AC_SSP5.2020_q1) / shape_ely_diff_m$cons_AC_SSP5.2020_q1
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1) | shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 = ifelse(shape_ely_diff_m$cons_AC_SSP5.2050_q1==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1)

shape_ely_diff_m$cons_AC_SSP5.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP5.2020_q3)
shape_ely_diff_m$cons_AC_SSP5.2050_q3 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2050_q3), 0, shape_ely_diff_m$cons_AC_SSP5.2050_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 = (shape_ely_diff_m$cons_AC_SSP5.2050_q3 - shape_ely_diff_m$cons_AC_SSP5.2020_q3) / shape_ely_diff_m$cons_AC_SSP5.2020_q3
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3) | shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 = ifelse(shape_ely_diff_m$cons_AC_SSP5.2050_q3==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3)

shape_ely_diff_m$cons_AC_SSP1.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020), 0, shape_ely_diff_m$cons_AC_SSP1.2020)
shape_ely_diff_m$cons_AC_SSP1.2050 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2050), 0, shape_ely_diff_m$cons_AC_SSP1.2050)
shape_ely_diff_m$pct_chg_cons_AC_SSP1 = (shape_ely_diff_m$cons_AC_SSP1.2050 - shape_ely_diff_m$cons_AC_SSP1.2020) / shape_ely_diff_m$cons_AC_SSP1.2020
shape_ely_diff_m$pct_chg_cons_AC_SSP1 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP1) | shape_ely_diff_m$pct_chg_cons_AC_SSP1 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1 = ifelse(shape_ely_diff_m$cons_AC_SSP1.2050==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP1<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP1)

shape_ely_diff_m$cons_AC_SSP1.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP1.2020_q1)
shape_ely_diff_m$cons_AC_SSP1.2050_q1 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2050_q1), 0, shape_ely_diff_m$cons_AC_SSP1.2050_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 = (shape_ely_diff_m$cons_AC_SSP1.2050_q1 - shape_ely_diff_m$cons_AC_SSP1.2020_q1) / shape_ely_diff_m$cons_AC_SSP1.2020_q1
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1) | shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 = ifelse(shape_ely_diff_m$cons_AC_SSP1.2050_q1==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1)

shape_ely_diff_m$cons_AC_SSP1.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP1.2020_q3)
shape_ely_diff_m$cons_AC_SSP1.2050_q3 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2050_q3), 0, shape_ely_diff_m$cons_AC_SSP1.2050_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 = (shape_ely_diff_m$cons_AC_SSP1.2050_q3 - shape_ely_diff_m$cons_AC_SSP1.2020_q3) / shape_ely_diff_m$cons_AC_SSP1.2020_q3
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3) | shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 = ifelse(shape_ely_diff_m$cons_AC_SSP1.2050_q3==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3)

shape_ely_diff_m$cons_AC_SSP3.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020), 0, shape_ely_diff_m$cons_AC_SSP3.2020)
shape_ely_diff_m$cons_AC_SSP3.2050 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2050), 0, shape_ely_diff_m$cons_AC_SSP3.2050)
shape_ely_diff_m$pct_chg_cons_AC_SSP3 = (shape_ely_diff_m$cons_AC_SSP3.2050 - shape_ely_diff_m$cons_AC_SSP3.2020) / shape_ely_diff_m$cons_AC_SSP3.2020
shape_ely_diff_m$pct_chg_cons_AC_SSP3 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP3) | shape_ely_diff_m$pct_chg_cons_AC_SSP3 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3 = ifelse(shape_ely_diff_m$cons_AC_SSP3.2050==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP3<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP3)

shape_ely_diff_m$cons_AC_SSP3.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP3.2020_q1)
shape_ely_diff_m$cons_AC_SSP3.2050_q1 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2050_q1), 0, shape_ely_diff_m$cons_AC_SSP3.2050_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 = (shape_ely_diff_m$cons_AC_SSP3.2050_q1 - shape_ely_diff_m$cons_AC_SSP3.2020_q1) / shape_ely_diff_m$cons_AC_SSP3.2020_q1
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1) | shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 = ifelse(shape_ely_diff_m$cons_AC_SSP3.2050_q1==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1)

shape_ely_diff_m$cons_AC_SSP3.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP3.2020_q3)
shape_ely_diff_m$cons_AC_SSP3.2050_q3 =ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2050_q3), 0, shape_ely_diff_m$cons_AC_SSP3.2050_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 = (shape_ely_diff_m$cons_AC_SSP3.2050_q3 - shape_ely_diff_m$cons_AC_SSP3.2020_q3) / shape_ely_diff_m$cons_AC_SSP3.2020_q3
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 = ifelse(!is.finite(shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3) | shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 >1, 1, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 = ifelse(shape_ely_diff_m$cons_AC_SSP3.2050_q3==0, 0, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 = ifelse(shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3<=0, NA, shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3)

####

drg1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP2, pct_chg_cons_AC_SSP2_q1, pct_chg_cons_AC_SSP2_q3, pop_ssps_data_ssp2_2050, SSP2.2050, id)
drg1 = arrange(drg1, pct_chg_cons_AC_SSP2)

drg2 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP5, pct_chg_cons_AC_SSP5_q1, pct_chg_cons_AC_SSP5_q3, pop_ssps_data_ssp5_2050, SSP5.2050, id)
drg2 = arrange(drg2, pct_chg_cons_AC_SSP5)

drg1 <- filter(drg1, pct_chg_cons_AC_SSP2>0)
drg2 <- filter(drg2, pct_chg_cons_AC_SSP5>0)

drg3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP1, pct_chg_cons_AC_SSP1_q1, pct_chg_cons_AC_SSP1_q3, pop_ssps_data_ssp1_2050, SSP1.2050, id)
drg3 = arrange(drg3, pct_chg_cons_AC_SSP1)

drg4 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP3, pct_chg_cons_AC_SSP3_q1, pct_chg_cons_AC_SSP3_q3, pop_ssps_data_ssp3_2050, SSP3.2050, id)
drg4 = arrange(drg4, pct_chg_cons_AC_SSP3)

drg3 <- filter(drg3, pct_chg_cons_AC_SSP1>0)
drg4 <- filter(drg4, pct_chg_cons_AC_SSP3>0)

####

a = ggplot()+
  theme_classic()+
  geom_line(data=drg3, aes(x=pct_chg_cons_AC_SSP1*100, y=(cumsum(pop_ssps_data_ssp1_2050*SSP1.2050)/1e9), colour="SSP126"))+
  geom_line(data=drg1, aes(x=pct_chg_cons_AC_SSP2*100, y=(cumsum(pop_ssps_data_ssp2_2050*SSP2.2050)/1e9), colour="SSP245"))+
  geom_line(data=drg4, aes(x=pct_chg_cons_AC_SSP3*100, y=(cumsum(pop_ssps_data_ssp3_2050*SSP3.2050)/1e9), colour="SSP370"))+
  geom_line(data=drg2, aes(x=pct_chg_cons_AC_SSP5*100, y=(cumsum(pop_ssps_data_ssp5_2050*SSP5.2050)/1e9), colour="SSP585"))+
  xlab("Percentage (%) change in AC-driven \nelectricity consumption")+
  ylab("People with AC (cumulative, billions)")+
  scale_colour_manual(name="", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  scale_fill_manual(name="", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), labels=c("0", "25", "50", "75", "100+"))+
  theme(legend.position = "none", legend.direction = "horizontal")

###

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

##
# calculate IQR of GCMs

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ac <- bind_cols(shape_ac, future_acc_q1, future_acc_q3)

#####

shape_ely_diff <- shape_ely_diff %>% dplyr::select(id, contains(".20") & contains("SSP"))
shape_ely_diff$region <- NULL
colnames(shape_ely_diff)[2:261] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:261])

shape_ely_diff <- ungroup(shape_ely_diff)

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ely_diff <- bind_cols(shape_ely_diff, future_acc_q1, future_acc_q3)

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff, "id")
shape_ely_diff_m[is.na(shape_ely_diff_m)] = 0

# x axis: percentage change in total energy use (from negative to >100)
# y axis: million people (cumulative sum)
# facets: region
# lines coloured by SSP/decile

shape_ely_diff_m$cons_AC_SSP2.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020), 0, shape_ely_diff_m$cons_AC_SSP2.2020)
shape_ely_diff_m$pct_chg_cons_AC_SSP2 = (shape_ely_diff_m$cons_AC_SSP2.2050 - shape_ely_diff_m$cons_AC_SSP2.2020) 

shape_ely_diff_m$cons_AC_SSP5.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020), 0, shape_ely_diff_m$cons_AC_SSP5.2020)
shape_ely_diff_m$pct_chg_cons_AC_SSP5 = (shape_ely_diff_m$cons_AC_SSP5.2050 - shape_ely_diff_m$cons_AC_SSP5.2020) 

shape_ely_diff_m$cons_AC_SSP1.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020), 0, shape_ely_diff_m$cons_AC_SSP1.2020)
shape_ely_diff_m$pct_chg_cons_AC_SSP1 = (shape_ely_diff_m$cons_AC_SSP1.2050 - shape_ely_diff_m$cons_AC_SSP1.2020) 

shape_ely_diff_m$cons_AC_SSP3.2020 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020), 0, shape_ely_diff_m$cons_AC_SSP3.2020)
shape_ely_diff_m$pct_chg_cons_AC_SSP3 = (shape_ely_diff_m$cons_AC_SSP3.2050 - shape_ely_diff_m$cons_AC_SSP3.2020) 

drg1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP2, pop_ssps_data_ssp2_2050, id)
drg1 <- merge(drg1, shape_ac %>% dplyr::select(id, SSP2.2050), "id")

drg2 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP5, pop_ssps_data_ssp5_2050, id)
drg2 <- merge(drg2, shape_ac%>% dplyr::select(id, SSP5.2050), "id")

drg1 <- filter(drg1, pct_chg_cons_AC_SSP2>0)
drg2 <- filter(drg2, pct_chg_cons_AC_SSP5>0)

drg1 = arrange(drg1, pct_chg_cons_AC_SSP2)
drg2 = arrange(drg2, pct_chg_cons_AC_SSP5)

drg1$cs <- cumsum(replace_na(drg1$pop_ssps_data_ssp2_2050*drg1$SSP2.2050, 0))/1e9
drg2$cs <- cumsum(replace_na(drg2$pop_ssps_data_ssp5_2050*drg2$SSP5.2050, 0))/1e9

drg3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP1, pop_ssps_data_ssp1_2050, id)
drg3 <- merge(drg3, shape_ac %>% dplyr::select(id, SSP1.2050), "id")

drg4 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP3, pop_ssps_data_ssp3_2050, id)
drg4 <- merge(drg4, shape_ac%>% dplyr::select(id, SSP3.2050), "id")

drg3 <- filter(drg3, pct_chg_cons_AC_SSP1>0)
drg4 <- filter(drg4, pct_chg_cons_AC_SSP3>0)

drg3 = arrange(drg3, pct_chg_cons_AC_SSP1)
drg4 = arrange(drg4, pct_chg_cons_AC_SSP3)

drg3$cs <- cumsum(replace_na(drg3$pop_ssps_data_ssp1_2050*drg3$SSP1.2050, 0))/1e9
drg4$cs <- cumsum(replace_na(drg4$pop_ssps_data_ssp3_2050*drg4$SSP3.2050, 0))/1e9

############

shape_ely_diff_m$cons_AC_SSP2.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP2.2020_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q1 = (shape_ely_diff_m$cons_AC_SSP2.2050_q1 - shape_ely_diff_m$cons_AC_SSP2.2020_q1) 

shape_ely_diff_m$cons_AC_SSP5.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP5.2020_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q1 = (shape_ely_diff_m$cons_AC_SSP5.2050_q1 - shape_ely_diff_m$cons_AC_SSP5.2020_q1) 

drg1_q1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP2_q1, pop_ssps_data_ssp2_2050, id)
drg1_q1 <- merge(drg1_q1, shape_ac %>% dplyr::select(id, SSP2.2050_q1), "id")

drg2_q1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP5_q1, pop_ssps_data_ssp5_2050, id)
drg2_q1 <- merge(drg2_q1, shape_ac %>% dplyr::select(id, SSP5.2050_q1), "id")

drg1_q1 <- filter(drg1_q1, pct_chg_cons_AC_SSP2_q1>0)
drg2_q1 <- filter(drg2_q1, pct_chg_cons_AC_SSP5_q1>0)

drg1_q1 = arrange(drg1_q1, pct_chg_cons_AC_SSP2_q1)
drg2_q1 = arrange(drg2_q1, pct_chg_cons_AC_SSP5_q1)

drg1_q1$cs_q1 <- cumsum(replace_na(drg1_q1$pop_ssps_data_ssp2_2050*drg1_q1$SSP2.2050_q1, 0))/1e9
drg2_q1$cs_q1 <- cumsum(replace_na(drg2_q1$pop_ssps_data_ssp5_2050*drg2_q1$SSP5.2050_q1, 0))/1e9

######

shape_ely_diff_m$cons_AC_SSP2.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP2.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP2.2020_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP2_q3 = (shape_ely_diff_m$cons_AC_SSP2.2050_q3 - shape_ely_diff_m$cons_AC_SSP2.2020_q3) 

shape_ely_diff_m$cons_AC_SSP5.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP5.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP5.2020_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP5_q3 = (shape_ely_diff_m$cons_AC_SSP5.2050_q3 - shape_ely_diff_m$cons_AC_SSP5.2020_q3) 

drg1_q3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP2_q3, pop_ssps_data_ssp2_2050, id)
drg1_q3 <- merge(drg1_q3, shape_ac %>% dplyr::select(id, SSP2.2050_q3), "id")

drg2_q3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP5_q3, pop_ssps_data_ssp5_2050, id)
drg2_q3 <- merge(drg2_q3, shape_ac %>% dplyr::select(id, SSP5.2050_q3), "id")

drg1_q3 <- filter(drg1_q3, pct_chg_cons_AC_SSP2_q3>0)
drg2_q3 <- filter(drg2_q3, pct_chg_cons_AC_SSP5_q3>0)

drg1_q3 = arrange(drg1_q3, pct_chg_cons_AC_SSP2_q3)
drg2_q3 = arrange(drg2_q3, pct_chg_cons_AC_SSP5_q3)

drg1_q3$cs_q3 <- cumsum(replace_na(drg1_q3$pop_ssps_data_ssp2_2050*drg1_q3$SSP2.2050_q3, 0))/1e9
drg2_q3$cs_q3 <- cumsum(replace_na(drg2_q3$pop_ssps_data_ssp5_2050*drg2_q3$SSP5.2050_q3, 0))/1e9

####

drg1_qq <- merge(drg1_q1, drg1_q3, "id", all=T)
drg2_qq <- merge(drg2_q1, drg2_q3, "id", all=T)

drg1_qq$pct_chg_cons_AC_SSP2_q1 <- (drg1_qq$pct_chg_cons_AC_SSP2_q1 + drg1_qq$pct_chg_cons_AC_SSP2_q3)/2

drg2_qq$pct_chg_cons_AC_SSP5_q1 <- (drg2_qq$pct_chg_cons_AC_SSP5_q1 + drg2_qq$pct_chg_cons_AC_SSP5_q3)/2


library(modelbased)

drg1_qq <- na.omit(drg1_qq)
drg2_qq <- na.omit(drg2_qq)

library(mgcv)

drg1_qq$cs_q1 <- predict(gam(drg1_qq$cs_q1 ~ s(drg1_qq$pct_chg_cons_AC_SSP2_q1, bs = "cs")))
drg1_qq$cs_q3 <- predict(gam(drg1_qq$cs_q3 ~ s(drg1_qq$pct_chg_cons_AC_SSP2_q1, bs = "cs")))

drg2_qq$cs_q1 <- predict(gam(drg2_qq$cs_q1 ~ s(drg2_qq$pct_chg_cons_AC_SSP5_q1, bs = "cs")))
drg2_qq$cs_q3 <- predict(gam(drg2_qq$cs_q3 ~ s(drg2_qq$pct_chg_cons_AC_SSP5_q1, bs = "cs")))

shape_ely_diff_m$cons_AC_SSP1.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP1.2020_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q1 = (shape_ely_diff_m$cons_AC_SSP1.2050_q1 - shape_ely_diff_m$cons_AC_SSP1.2020_q1) 

shape_ely_diff_m$cons_AC_SSP3.2020_q1 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020_q1), 0, shape_ely_diff_m$cons_AC_SSP3.2020_q1)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q1 = (shape_ely_diff_m$cons_AC_SSP3.2050_q1 - shape_ely_diff_m$cons_AC_SSP3.2020_q1) 

drg3_q1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP1_q1, pop_ssps_data_ssp1_2050, id)
drg3_q1 <- merge(drg3_q1, shape_ac %>% dplyr::select(id, SSP1.2050_q1), "id")

drg4_q1 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP3_q1, pop_ssps_data_ssp3_2050, id)
drg4_q1 <- merge(drg4_q1, shape_ac %>% dplyr::select(id, SSP3.2050_q1), "id")

drg3_q1 <- filter(drg3_q1, pct_chg_cons_AC_SSP1_q1>0)
drg4_q1 <- filter(drg4_q1, pct_chg_cons_AC_SSP3_q1>0)

drg3_q1 = arrange(drg3_q1, pct_chg_cons_AC_SSP1_q1)
drg4_q1 = arrange(drg4_q1, pct_chg_cons_AC_SSP3_q1)

drg3_q1$cs_q1 <- cumsum(replace_na(drg3_q1$pop_ssps_data_ssp1_2050*drg3_q1$SSP1.2050_q1, 0))/1e9
drg4_q1$cs_q1 <- cumsum(replace_na(drg4_q1$pop_ssps_data_ssp3_2050*drg4_q1$SSP3.2050_q1, 0))/1e9

######

shape_ely_diff_m$cons_AC_SSP1.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP1.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP1.2020_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP1_q3 = (shape_ely_diff_m$cons_AC_SSP1.2050_q3 - shape_ely_diff_m$cons_AC_SSP1.2020_q3) 

shape_ely_diff_m$cons_AC_SSP3.2020_q3 = ifelse(is.na(shape_ely_diff_m$cons_AC_SSP3.2020_q3), 0, shape_ely_diff_m$cons_AC_SSP3.2020_q3)
shape_ely_diff_m$pct_chg_cons_AC_SSP3_q3 = (shape_ely_diff_m$cons_AC_SSP3.2050_q3 - shape_ely_diff_m$cons_AC_SSP3.2020_q3) 

drg3_q3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP1_q3, pop_ssps_data_ssp1_2050, id)
drg3_q3 <- merge(drg3_q3, shape_ac %>% dplyr::select(id, SSP1.2050_q3), "id")

drg4_q3 = shape_ely_diff_m %>% dplyr::select(pct_chg_cons_AC_SSP3_q3, pop_ssps_data_ssp3_2050, id)
drg4_q3 <- merge(drg4_q3, shape_ac %>% dplyr::select(id, SSP3.2050_q3), "id")

drg3_q3 <- filter(drg3_q3, pct_chg_cons_AC_SSP1_q3>0)
drg4_q3 <- filter(drg4_q3, pct_chg_cons_AC_SSP3_q3>0)

drg3_q3 = arrange(drg3_q3, pct_chg_cons_AC_SSP1_q3)
drg4_q3 = arrange(drg4_q3, pct_chg_cons_AC_SSP3_q3)

drg3_q3$cs_q3 <- cumsum(replace_na(drg3_q3$pop_ssps_data_ssp1_2050*drg3_q3$SSP1.2050_q3, 0))/1e9
drg4_q3$cs_q3 <- cumsum(replace_na(drg4_q3$pop_ssps_data_ssp3_2050*drg4_q3$SSP3.2050_q3, 0))/1e9

####

drg3_qq <- merge(drg3_q1, drg3_q3, "id", all=T)
drg4_qq <- merge(drg4_q1, drg4_q3, "id", all=T)

drg3_qq$pct_chg_cons_AC_SSP1_q1 <- (drg3_qq$pct_chg_cons_AC_SSP1_q1 + drg3_qq$pct_chg_cons_AC_SSP1_q3)/2

drg4_qq$pct_chg_cons_AC_SSP3_q1 <- (drg4_qq$pct_chg_cons_AC_SSP3_q1 + drg4_qq$pct_chg_cons_AC_SSP3_q3)/2


library(modelbased)

drg3_qq <- na.omit(drg3_qq)
drg4_qq <- na.omit(drg4_qq)

library(mgcv)

drg3_qq$cs_q1 <- predict(gam(drg3_qq$cs_q1 ~ s(drg3_qq$pct_chg_cons_AC_SSP1_q1, bs = "cs")))
drg3_qq$cs_q3 <- predict(gam(drg3_qq$cs_q3 ~ s(drg3_qq$pct_chg_cons_AC_SSP1_q1, bs = "cs")))

drg4_qq$cs_q1 <- predict(gam(drg4_qq$cs_q1 ~ s(drg4_qq$pct_chg_cons_AC_SSP3_q1, bs = "cs")))
drg4_qq$cs_q3 <- predict(gam(drg4_qq$cs_q3 ~ s(drg4_qq$pct_chg_cons_AC_SSP3_q1, bs = "cs")))

b = ggplot()+
  theme_classic()+
  geom_line(data=drg3, aes(x=pct_chg_cons_AC_SSP1, y=cs, colour="SSP126"))+
  geom_ribbon(data=drg3_qq, aes(x=pct_chg_cons_AC_SSP1_q1, ymin=cs_q1, ymax=cs_q3, fill="SSP126"), alpha=0.1)+
  geom_line(data=drg1, aes(x=pct_chg_cons_AC_SSP2, y=cs, colour="SSP245"))+
  geom_ribbon(data=drg1_qq, aes(x=pct_chg_cons_AC_SSP2_q1, ymin=cs_q1, ymax=cs_q3, fill="SSP245"), alpha=0.1)+
  geom_line(data=drg4, aes(x=pct_chg_cons_AC_SSP3, y=cs, colour="SSP370"))+
  geom_ribbon(data=drg4_qq, aes(x=pct_chg_cons_AC_SSP3_q1, ymin=cs_q1, ymax=cs_q3, fill="SSP370"), alpha=0.1)+
  geom_line(data=drg2, aes(x=pct_chg_cons_AC_SSP5, y=cs, colour="SSP585"))+
  geom_ribbon(data=drg2_qq, aes(x=pct_chg_cons_AC_SSP5_q1, ymin=cs_q1, ymax=cs_q3, fill="SSP585"), alpha=0.1)+
  xlab("Change in AC-driven electricity consumption (kWh/hh/yr), 2020-2050")+
  ylab("Population with AC (cumulative, billions)")+
  scale_colour_manual(name="", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  scale_fill_manual(name="", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  xlim(0, 3000)+
  theme(legend.position = "bottom", legend.direction = "horizontal")

library(patchwork)

(a + b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave("results/graphs_tables/cum_line_abs_ac.pdf", b, scale=1.2, width = 7, height = 4.5)
  

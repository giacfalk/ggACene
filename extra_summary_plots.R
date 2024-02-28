
#######

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

## hhsize

shape_ac[,grep("pop_", colnames(shape_ac))] <- shape_ac[,grep("pop_", colnames(shape_ac))] / shape_ac$hhsize

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

##########

shape_ac_s <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2010_q1=weighted.mean(SSP2.2010_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q1=weighted.mean(SSP2.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q1=weighted.mean(SSP5.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2010_q3=weighted.mean(SSP2.2010_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q3=weighted.mean(SSP2.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q3=weighted.mean(SSP5.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2010_q1=weighted.mean(SSP1.2010_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050_q1=weighted.mean(SSP1.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050_q1=weighted.mean(SSP3.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2010_q3=weighted.mean(SSP1.2010_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050_q3=weighted.mean(SSP1.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050_q3=weighted.mean(SSP3.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T))

library(acid)

shape_ac_s_gini <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::reframe(SSP2.2010=weighted.gini(SSP2.2010, pop_ssps_data_ssp2_2020), SSP2.2050=weighted.gini(SSP2.2050, pop_ssps_data_ssp2_2050), SSP5.2050=weighted.gini(SSP5.2050, pop_ssps_data_ssp2_2050), SSP2.2010_q1=weighted.gini(SSP2.2010_q1, pop_ssps_data_ssp2_2020), SSP2.2050_q1=weighted.gini(SSP2.2050_q1, pop_ssps_data_ssp2_2050), SSP5.2050_q1=weighted.gini(SSP5.2050_q1, pop_ssps_data_ssp2_2050), SSP2.2010_q3=weighted.gini(SSP2.2010_q3, pop_ssps_data_ssp2_2020), SSP2.2050_q3=weighted.gini(SSP2.2050_q3, pop_ssps_data_ssp2_2050), SSP5.2050_q3=weighted.gini(SSP5.2050_q3, pop_ssps_data_ssp2_2050), SSP1.2010=weighted.gini(SSP1.2010, pop_ssps_data_ssp2_2020), SSP1.2050=weighted.gini(SSP1.2050, pop_ssps_data_ssp2_2050), SSP3.2050=weighted.gini(SSP3.2050, pop_ssps_data_ssp2_2050), SSP1.2010_q1=weighted.gini(SSP1.2010_q1, pop_ssps_data_ssp2_2020), SSP1.2050_q1=weighted.gini(SSP1.2050_q1, pop_ssps_data_ssp2_2050), SSP3.2050_q1=weighted.gini(SSP3.2050_q1, pop_ssps_data_ssp2_2050), SSP1.2010_q3=weighted.gini(SSP1.2010_q3, pop_ssps_data_ssp2_2020), SSP1.2050_q3=weighted.gini(SSP1.2050_q3, pop_ssps_data_ssp2_2050), SSP3.2050_q3=weighted.gini(SSP3.2050_q3, pop_ssps_data_ssp2_2050))

cols.num <- colnames(shape_ac_s_gini)[2:19]

shape_ac_s_gini[cols.num] <- sapply(shape_ac_s_gini[cols.num],as.numeric)

shape_ac_s_gini <- shape_ac_s_gini %>% group_by(ISO3) %>% dplyr::summarise_at( .vars = colnames(.)[2:19], mean, na.rm=T)

colnames(shape_ac_s_gini)[2:19] <- paste0(colnames(shape_ac_s_gini)[2:19], "_gini")

#####################

shape_ely_diff <- shape_ely_diff %>% dplyr::select(contains(".20") & contains("SSP"), id)
shape_ely_diff$ISO3 <- NULL
colnames(shape_ely_diff)[2:261] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:261])

shape_ely_diff <- ungroup(shape_ely_diff)

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y)) & (!grepl("GDP", names(shape_ely_diff)) & !grepl("URB", names(shape_ely_diff)) & !grepl("women", names(shape_ely_diff)) & !grepl("age", names(shape_ely_diff)) & !grepl("edu", names(shape_ely_diff)))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y)) & (!grepl("GDP", names(shape_ely_diff)) & !grepl("URB", names(shape_ely_diff)) & !grepl("women", names(shape_ely_diff)) & !grepl("age", names(shape_ely_diff)) & !grepl("edu", names(shape_ely_diff)))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ely_diff <- bind_cols(shape_ely_diff, future_acc_q1, future_acc_q3)

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff, "id")

##########

shape_ely_diff_s <- dplyr::group_by(shape_ely_diff_m, ISO3) %>% dplyr::summarise(ely_total_SSP2_2010= sum(cons_AC_SSP2.2010*(SSP2.2010*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050= sum(cons_AC_SSP2.2050*(SSP2.2050*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050= sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2010_q1= sum(cons_AC_SSP2.2010_q1*(SSP2.2010_q1*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q1= sum(cons_AC_SSP2.2050_q1*(SSP2.2050_q1*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q1= sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2010_q3= sum(cons_AC_SSP2.2010_q3*(SSP2.2010_q3*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q3= sum(cons_AC_SSP2.2050_q3*(SSP2.2050_q3*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q3= sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP1_2010= sum(cons_AC_SSP1.2010*(SSP1.2010*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050= sum(cons_AC_SSP1.2050*(SSP1.2050*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050= sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9, ely_total_SSP1_2010_q1= sum(cons_AC_SSP1.2010_q1*(SSP1.2010_q1*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050_q1= sum(cons_AC_SSP1.2050_q1*(SSP1.2050_q1*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050_q1= sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9, ely_total_SSP1_2010_q3= sum(cons_AC_SSP1.2010_q3*(SSP1.2010_q3*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050_q3= sum(cons_AC_SSP1.2050_q3*(SSP1.2050_q3*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050_q3= sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9)

#

shape_ely_diff_gini <- dplyr::group_by(na.omit(shape_ely_diff_m), ISO3) %>% dplyr::reframe(ely_total_SSP2_2010= weighted.gini(cons_AC_SSP2.2010, (SSP2.2010*(pop_ssps_data_ssp2_2020))), ely_total_SSP2_2050= weighted.gini(cons_AC_SSP2.2050, (SSP2.2050*(pop_ssps_data_ssp2_2050))), ely_total_SSP5_2050= weighted.gini(cons_AC_SSP5.2050, (SSP5.2050*(pop_ssps_data_ssp5_2050))), ely_total_SSP1_2050= weighted.gini(cons_AC_SSP1.2050, (SSP1.2050*(pop_ssps_data_ssp1_2050))), ely_total_SSP3_2050= weighted.gini(cons_AC_SSP3.2050, (SSP3.2050*(pop_ssps_data_ssp3_2050))))
#

cols.num <- colnames(shape_ely_diff_gini)[2:6]

shape_ely_diff_gini[cols.num] <- sapply(shape_ely_diff_gini[cols.num],as.numeric)

shape_ely_diff_gini <- shape_ely_diff_gini %>% group_by(ISO3) %>% dplyr::summarise_at( .vars = colnames(.)[2:6], mean, na.rm=T)

colnames(shape_ely_diff_gini)[2:6] <- paste0(colnames(shape_ely_diff_gini)[2:6], "_gini")

##################################
##################################
##################################

acglobtot = shape_ac %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T))

elyglobtot =  dplyr::group_by(shape_ely_diff_m) %>% dplyr::summarise(cons_AC_SSP2.2010=weighted.mean(cons_AC_SSP2.2010, (pop_ssps_data_ssp2_2020*SSP2.2010), na.rm=T), cons_AC_SSP2.2050=weighted.mean(cons_AC_SSP2.2050, (pop_ssps_data_ssp2_2050*SSP2.2050), na.rm=T))

acglobtot$SSP2.2050 <- acglobtot$SSP2.2050 / acglobtot$SSP2.2010
elyglobtot$cons_AC_SSP2.2050 <- elyglobtot$cons_AC_SSP2.2050 / elyglobtot$cons_AC_SSP2.2010

####

shape_ely_diff_s <- dplyr::group_by(shape_ely_diff_m, ISO3) %>% dplyr::summarise(cons_AC_SSP2.2010=weighted.mean(cons_AC_SSP2.2010, pop_ssps_data_ssp2_2020*SSP2.2010, na.rm=T), cons_AC_SSP2.2050=weighted.mean(cons_AC_SSP2.2050, pop_ssps_data_ssp2_2050*SSP2.2050, na.rm=T))

ss1 <- dplyr::select(shape_ac_s, 1, 3)
ss1$geometry <- NULL
ss1$SSP2.2050 <- ss1$SSP2.2050 / shape_ac_s$SSP2.2010 

ss1$SSP2.2050 <- ss1$SSP2.2050 / acglobtot$SSP2.2050
  
ss2 <- dplyr::select(shape_ely_diff_s, 3)
ss2$geometry <- NULL
ss2$cons_AC_SSP2.2050 <- ss2$cons_AC_SSP2.2050 / shape_ely_diff_s$cons_AC_SSP2.2010

ss2$cons_AC_SSP2.2050 <- ss2$cons_AC_SSP2.2050 / elyglobtot$cons_AC_SSP2.2050

sss <- bind_cols(ss1, ss2)
sss$region <- countrycode::countrycode(sss$ISO3, 'iso3c', 'region')

###

ziocountries <- tail(arrange(wrld_simpl, POP2005), 30)$ISO3

library(ggrepel)

ggplot(sss)+
  theme_classic()+
  geom_vline(xintercept = 1, alpha=0.5)+
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_point(aes(x=SSP2.2050, y=cons_AC_SSP2.2050, colour=region))+
  xlab("Relative change in national AC penetration rate (%) compared to global avg. change")+
  ylab("Relative change in avg. AC electricity consumption (kWh/hh/yr.) compared to global avg. change")+
  ggtitle("Delta to 2050, SSP245")+
  scale_colour_brewer(palette = "Set1")+
  scale_x_continuous(limits=c(0,3))+
  scale_y_continuous(limits=c(0,3))+
  geom_text(aes(x=2.5, y = 3), label="Faster growth in extensive \nand intensive margins")+
  geom_text(aes(x=0.4, y = 3), label="Slower growth in extensive margin \nbut faster growth in intensive margin")+
  geom_text(aes(x=2.5, y = 0.3), label="Faster growth in extensive margin \n but slower growth in intensive margin")+
  geom_text(aes(x=0.35, y = 0.3), label="Slower growth in extensive \nand intensive margins")+
  geom_label_repel(data = sss %>% 
                    mutate(label = ifelse(ISO3 %in% ziocountries, as.character(ISO3), "")),
                  aes(x=SSP2.2050, y=cons_AC_SSP2.2050, label = label), 
                  box.padding = 1,
                  show.legend = FALSE, max.overlaps =Inf) 

ggsave("results/graphs_tables/scatter_changes_ssp2.png", scale=2.35, height=3, width = 4.5)

##################################
##################################
##################################


#faceted by region, ac penetration by decile of income (today and in two scenarios)

shape_ac <- shape_ac %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

shape_ac_d1 <- filter(shape_ac, decile==1)
shape_ac_d5 <- filter(shape_ac, decile==5)

ggplot()+
  theme_classic()+
  geom_histogram(data=shape_ac_d1 %>% filter((SSP2.2050-SSP2.2010)>0), aes(x=(SSP2.2050-SSP2.2010)*100, weight=pop_ssps_data_ssp2_2050, fill="Q1"), alpha=0.5)+
  geom_histogram(data=shape_ac_d5 %>% filter((SSP2.2050-SSP2.2010)>0), aes(x=(SSP2.2050-SSP2.2010)*100, weight=pop_ssps_data_ssp2_2050, fill="Q5"), alpha=0.5)+
  facet_wrap(vars(region), scales = "free")+
  ggtitle("")+
  xlab("Absolute change in the AC penetration rate")+
  ylab("")+
  scale_fill_discrete(name="Income quintile")

ggsave("results/graphs_tables/hist_quintiles_ac_ssp2.png", scale=2.35, height=3, width = 4.5)

###

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff, "id")
shape_ely_diff_m <- shape_ely_diff_m %>% group_by(region.x) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

shape_ely_d1 <- filter(shape_ely_diff_m, decile==1)
shape_ely_d5 <- filter(shape_ely_diff_m, decile==5)

ggplot()+
  theme_classic()+
  geom_histogram(data=shape_ely_d1 %>% filter((cons_AC_SSP2.2050-cons_AC_SSP2.2010)>0), aes(x=(cons_AC_SSP2.2050-cons_AC_SSP2.2010), weight=pop_ssps_data_ssp2_2050*SSP2.2050, fill="Q1"), alpha=0.5)+
  geom_histogram(data=shape_ely_d5 %>% filter((cons_AC_SSP2.2050-cons_AC_SSP2.2010)>0), aes(x=(cons_AC_SSP2.2050-cons_AC_SSP2.2010), weight=pop_ssps_data_ssp2_2050, fill="Q5"), alpha=0.5)+
  facet_wrap(vars(region.x), scales = "free")+
  ggtitle("")+
  xlab("Absolute change in AC electricity (kWh/hh/yr.)")+
  ylab("")+
  scale_fill_discrete(name="Income quintile")

ggsave("results/graphs_tables/hist_quintiles_ely_ssp2.png", scale=2.35, height=3, width = 4.5)

  
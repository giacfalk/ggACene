# stats for paper
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

acglobtot = dplyr::group_by(shape_ac) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

acglobtot

median(as.numeric(acglobtot)[3:6])

#


acglobtot = dplyr::group_by(shape_ac, region) %>% dplyr::summarise(SSP2.2010=sum(SSP2.2010* pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020=sum(SSP2.2020* pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=sum(SSP2.2050* pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=sum(SSP5.2050* pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2050=sum(SSP1.2050* pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2050=sum(SSP3.2050* pop_ssps_data_ssp3_2050, na.rm=T))

acglobtot <- bind_cols(acglobtot[,1], acglobtot[,c(2:7)] / 1e6)

acglobtot[,8] <- (rowMeans(acglobtot[,c(4:7)]) / acglobtot[,3]) - 1

colnames(acglobtot)[8] <- "growth"

arrange(acglobtot, growth)


acglobtot[,9] <- rowMeans(acglobtot[,c(4:7)]) - acglobtot[,3] 

colnames(acglobtot)[9] <- "abs"

arrange(acglobtot, abs)

dplyr::group_by(shape_ac, region) %>% dplyr::summarise(SSP1.2050=sum(pop_ssps_data_ssp1_2050, na.rm=T)/1e6)

#

acglobtotHH = dplyr::group_by(shape_ac) %>% dplyr::summarise(SSP2.2010=sum(SSP2.2010* pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020=sum(SSP2.2020* pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=sum(SSP2.2050* pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=sum(SSP5.2050* pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2050=sum(SSP1.2050* pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2050=sum(SSP3.2050* pop_ssps_data_ssp3_2050, na.rm=T))

acglobtotHH/1e9

median(as.numeric(acglobtotHH)[3:6])



#

shape_ely_diff <- shape_ely_diff %>% dplyr::select(starts_with("SSP"), id)
shape_ely_diff$ISO3 <- NULL
colnames(shape_ely_diff)[2:21] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:21])

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff, "id")

elyglobtot = dplyr::group_by(shape_ely_diff_m) %>% dplyr::summarise(ely_total_SSP2_2010= sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_ssp2_2020), na.rm = T)/1e9, ely_total_SSP2_2020= sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm = T)/1e9, ely_total_SSP2_2050= sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm = T)/1e9, ely_total_SSP5_2050= sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP1_2050= sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm = T)/1e9, ely_total_SSP3_2050= sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9)

elyglobtot

median(as.numeric(elyglobtot)[3:6])


#######

shape_ac <- shape_ac %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2020), 5))

paths <- dplyr::group_by(shape_ac, region, decile) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T), SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

paths_glob <- dplyr::group_by(shape_ac, decile) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T), SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

paths$type ="Region"
paths_glob$type ="Global"

paths <- bind_rows(paths, paths_glob)

paths <- reshape2::melt(paths, c(1, 2, 23))
paths$ssp <- substr(paths$variable, 1, 4)
paths$year <- as.numeric(substr(paths$variable, 6, 9))

paths$ISO3 <- as.character(paths$region)
paths$ISO3 <- ifelse(is.na(paths$ISO3), " GLOBAL", paths$ISO3)

paths$decile <- as.factor(paths$decile)

paths <- filter(paths, year==2020 | year==2050)

paths$level <- ifelse(paths$ssp=="SSP2" & paths$year==2020, "2020", ifelse(paths$ssp=="SSP2" & paths$year==2050, "2050, SSP245",  ifelse(paths$ssp=="SSP5" & paths$year==2050, "2050, SSP585", ifelse(paths$ssp=="SSP1" & paths$year==2050, "2050, SSP126", ifelse(paths$ssp=="SSP3" & paths$year==2050, "2050, SSP370", NA)))))

paths <- filter(paths, !is.na(level))
paths$level <- as.factor(paths$level)

paths %>% filter(type=="Global") %>% filter(decile==1)
paths %>% filter(type=="Global") %>% filter(decile==5)


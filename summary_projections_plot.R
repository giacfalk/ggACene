
#######

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
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

setwd(wd)

# descriptive statistics
#######

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

## hhsize

shape_ac[,grep("pop_", colnames(shape_ac))] <- shape_ac[,grep("pop_", colnames(shape_ac))] / shape_ac$hhsize

##

# some salient results

base <- (weighted.mean(shape_ely_diff$SSP2.2020, shape_ac$SSP2.2020*shape_ac$pop_ssps_data_ssp2_2020, na.rm=T)*((sum(shape_ac$pop_ssps_data_ssp2_2020*shape_ac$SSP2.2020))))/1e9
a <- (weighted.mean(shape_ely_diff$SSP2.2050, shape_ac$SSP2.2050*shape_ac$pop_ssps_data_ssp2_2050, na.rm=T)*((sum(shape_ac$pop_ssps_data_ssp2_2050*shape_ac$SSP2.2050))))/1e9
b <- (weighted.mean(shape_ely_diff$SSP5.2050, shape_ac$SSP5.2050*shape_ac$pop_ssps_data_ssp5_2050, na.rm=T)*((sum(shape_ac$pop_ssps_data_ssp5_2050*shape_ac$SSP5.2050))))/1e9
c <- (weighted.mean(shape_ely_diff$SSP1.2050, shape_ac$SSP1.2050*shape_ac$pop_ssps_data_ssp1_2050, na.rm=T)*((sum(shape_ac$pop_ssps_data_ssp1_2050*shape_ac$SSP1.2050))))/1e9
d <- (weighted.mean(shape_ely_diff$SSP3.2050, shape_ac$SSP3.2050*shape_ac$pop_ssps_data_ssp3_2050, na.rm=T)*((sum(shape_ac$pop_ssps_data_ssp3_2050*shape_ac$SSP3.2050))))/1e9

base
a
b
c
d

median(c(a, b, c, d))

#

weighted.mean(shape_ely_diff$SSP2.2020, shape_ac$SSP2.2020*shape_ac$pop_ssps_data_ssp2_2020, na.rm=T)
weighted.mean(shape_ely_diff$SSP2.2050, shape_ac$SSP2.2050*shape_ac$pop_ssps_data_ssp2_2050, na.rm=T)
weighted.mean(shape_ely_diff$SSP5.2050, shape_ac$SSP5.2050*shape_ac$pop_ssps_data_ssp5_2050, na.rm=T)
weighted.mean(shape_ely_diff$SSP1.2050, shape_ac$SSP1.2050*shape_ac$pop_ssps_data_ssp1_2050, na.rm=T)
weighted.mean(shape_ely_diff$SSP3.2050, shape_ac$SSP3.2050*shape_ac$pop_ssps_data_ssp3_2050, na.rm=T)

#

weighted.mean(shape_ely_diff$SSP2.2020[shape_ely_diff$region=="Sub-Saharan Africa"], shape_ac$SSP2.2020[shape_ely_diff$region=="Sub-Saharan Africa"]*shape_ac$pop_ssps_data_ssp2_2020[shape_ely_diff$region=="Sub-Saharan Africa"], na.rm=T)
weighted.mean(shape_ely_diff$SSP2.2020[shape_ely_diff$region=="North America"], shape_ac$SSP2.2020[shape_ely_diff$region=="North America"]*shape_ac$pop_ssps_data_ssp2_2020[shape_ely_diff$region=="North America"], na.rm=T)

#

base <- weighted.mean(shape_ac$SSP2.2020, shape_ac$pop_ssps_data_ssp2_2020, na.rm=T)
a <- weighted.mean(shape_ac$SSP2.2050, shape_ac$pop_ssps_data_ssp2_2050, na.rm=T)
b <- weighted.mean(shape_ac$SSP5.2050, shape_ac$pop_ssps_data_ssp5_2050, na.rm=T)
c <- weighted.mean(shape_ac$SSP1.2050, shape_ac$pop_ssps_data_ssp1_2050, na.rm=T)
d <- weighted.mean(shape_ac$SSP3.2050, shape_ac$pop_ssps_data_ssp3_2050, na.rm=T)

base * 100
a * 100
b * 100
c * 100
d * 100

median(c(a, b, c, d))


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

shape_ac_s <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2020_q1=weighted.mean(SSP2.2020_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q1=weighted.mean(SSP2.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q1=weighted.mean(SSP5.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2020_q3=weighted.mean(SSP2.2020_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q3=weighted.mean(SSP2.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q3=weighted.mean(SSP5.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2020_q1=weighted.mean(SSP1.2020_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050_q1=weighted.mean(SSP1.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050_q1=weighted.mean(SSP3.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP1.2020_q3=weighted.mean(SSP1.2020_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP1.2050_q3=weighted.mean(SSP1.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP3.2050_q3=weighted.mean(SSP3.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T))

library(acid)

shape_ac_s_gini <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::reframe(SSP2.2020=weighted.gini(SSP2.2020, pop_ssps_data_ssp2_2020), SSP2.2050=weighted.gini(SSP2.2050, pop_ssps_data_ssp2_2050), SSP5.2050=weighted.gini(SSP5.2050, pop_ssps_data_ssp2_2050), SSP2.2020_q1=weighted.gini(SSP2.2020_q1, pop_ssps_data_ssp2_2020), SSP2.2050_q1=weighted.gini(SSP2.2050_q1, pop_ssps_data_ssp2_2050), SSP5.2050_q1=weighted.gini(SSP5.2050_q1, pop_ssps_data_ssp2_2050), SSP2.2020_q3=weighted.gini(SSP2.2020_q3, pop_ssps_data_ssp2_2020), SSP2.2050_q3=weighted.gini(SSP2.2050_q3, pop_ssps_data_ssp2_2050), SSP5.2050_q3=weighted.gini(SSP5.2050_q3, pop_ssps_data_ssp2_2050), SSP1.2020=weighted.gini(SSP1.2020, pop_ssps_data_ssp2_2020), SSP1.2050=weighted.gini(SSP1.2050, pop_ssps_data_ssp2_2050), SSP3.2050=weighted.gini(SSP3.2050, pop_ssps_data_ssp2_2050), SSP1.2020_q1=weighted.gini(SSP1.2020_q1, pop_ssps_data_ssp2_2020), SSP1.2050_q1=weighted.gini(SSP1.2050_q1, pop_ssps_data_ssp2_2050), SSP3.2050_q1=weighted.gini(SSP3.2050_q1, pop_ssps_data_ssp2_2050), SSP1.2020_q3=weighted.gini(SSP1.2020_q3, pop_ssps_data_ssp2_2020), SSP1.2050_q3=weighted.gini(SSP1.2050_q3, pop_ssps_data_ssp2_2050), SSP3.2050_q3=weighted.gini(SSP3.2050_q3, pop_ssps_data_ssp2_2050))

cols.num <- colnames(shape_ac_s_gini)[2:19]

shape_ac_s_gini[cols.num] <- sapply(shape_ac_s_gini[cols.num],as.numeric)

shape_ac_s_gini <- shape_ac_s_gini %>% group_by(ISO3) %>% dplyr::summarise_at( .vars = colnames(.)[2:19], mean, na.rm=T)

colnames(shape_ac_s_gini)[2:19] <- paste0(colnames(shape_ac_s_gini)[2:19], "_gini")

#
library(maptools)
data("wrld_simpl")
wrld_simpl <- subset(wrld_simpl, wrld_simpl$ISO3!="ATA")
wrld_simpl <- st_as_sf(wrld_simpl)

shape_ac_s <- merge(shape_ac_s, wrld_simpl, "ISO3")
shape_ac_s <- st_as_sf(shape_ac_s)

shape_ac_s <- st_transform(shape_ac_s, "ESRI:54009")
wrld_simpl <- st_transform(wrld_simpl, "ESRI:54009")

# o2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ac_s, aes(fill=SSP2.2020*100), colour=NA, show.legend = F)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("National AC penetration rate, 2020")
# 
# 
# a2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ac_s, aes(fill=SSP2.2050*100), colour=NA, show.legend = F)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("National AC penetration rate, 2050, SSP2")
# 
# b2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ac_s, aes(fill=SSP5.2050*100), colour=NA)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("National AC penetration rate, 2050, SSP5")
# 
# library(patchwork)
# 
# a_p2_1 <- o2 + a2 + b2 + plot_annotation(tag_levels = list(c('A', 'B', 'C'))) +  plot_layout(guides = "collect", ncol=1) &   scale_fill_binned(name="%", type="viridis") & theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.width=unit(0.1,"npc"))

library(stargazer)

s_table <- shape_ac_s %>% dplyr::select(ISO3, SSP2.2020, SSP2.2050, SSP5.2050, SSP2.2020_q1, SSP2.2050_q1, SSP5.2050_q1, SSP2.2020_q3, SSP2.2050_q3, SSP5.2050_q3, SSP1.2020, SSP1.2050, SSP3.2050, SSP1.2020_q1, SSP1.2050_q1, SSP3.2050_q1, SSP1.2020_q3, SSP1.2050_q3, SSP3.2050_q3) %>% mutate(ISO3=as.character(ISO3))
s_table$geometry <- NULL

s_table <- bind_cols(s_table, shape_ac_s_gini %>% dplyr::select(-1))

s_table$SSP2.2020 <- paste0(round(s_table$SSP2.2020, 2), " (", round(s_table$SSP2.2020_q1, 2), " - ", round(s_table$SSP2.2020_q3, 2), ")")
s_table$SSP2.2050 <- paste0(round(s_table$SSP2.2050, 2), " (", round(s_table$SSP2.2050_q1, 2), " - ", round(s_table$SSP2.2050_q3, 2), ")")
s_table$SSP5.2050 <- paste0(round(s_table$SSP5.2050, 2), " (", round(s_table$SSP5.2050_q1, 2), " - ", round(s_table$SSP5.2050_q3, 2), ")")
s_table$SSP1.2050 <- paste0(round(s_table$SSP1.2050, 2), " (", round(s_table$SSP1.2050_q1, 2), " - ", round(s_table$SSP1.2050_q3, 2), ")")
s_table$SSP3.2050 <- paste0(round(s_table$SSP3.2050, 2), " (", round(s_table$SSP3.2050_q1, 2), " - ", round(s_table$SSP3.2050_q3, 2), ")")

s_table <- s_table %>% dplyr::select(-SSP2.2020_q1, -SSP2.2050_q1, -SSP5.2050_q1, -SSP2.2020_q3, -SSP2.2050_q3, -SSP5.2050_q3, -SSP1.2020_q1, -SSP1.2050_q1, -SSP3.2050_q1, -SSP1.2020_q3, -SSP1.2050_q3, -SSP3.2050_q3, -SSP2.2020_q1_gini, -SSP2.2050_q1_gini, -SSP5.2050_q1_gini, -SSP2.2020_q3_gini, -SSP2.2050_q3_gini, -SSP5.2050_q3_gini, -SSP1.2020_q1_gini, -SSP1.2050_q1_gini, -SSP3.2050_q1_gini, -SSP1.2020_q3_gini, -SSP1.2050_q3_gini, -SSP3.2050_q3_gini) %>% mutate(ISO3=as.character(ISO3))

s_table$SSP1.2020 <- NULL
s_table$SSP1.2020_gini <- NULL

s_table <- s_table[,c(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11)]

s_table <- s_table %>% dplyr::mutate_if(is.numeric, round, 2)

stargazer::stargazer(s_table, summary=F, out="results/graphs_tables/AC_penetrations.tex")

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

shape_ely_diff_s <- dplyr::group_by(shape_ely_diff_m, ISO3) %>% dplyr::summarise(ely_total_SSP2_2020= sum(cons_AC_SSP2.2020*(SSP2.2020*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050= sum(cons_AC_SSP2.2050*(SSP2.2050*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050= sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2020_q1= sum(cons_AC_SSP2.2020_q1*(SSP2.2020_q1*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q1= sum(cons_AC_SSP2.2050_q1*(SSP2.2050_q1*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q1= sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2020_q3= sum(cons_AC_SSP2.2020_q3*(SSP2.2020_q3*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q3= sum(cons_AC_SSP2.2050_q3*(SSP2.2050_q3*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q3= sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP1_2020= sum(cons_AC_SSP1.2020*(SSP1.2020*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050= sum(cons_AC_SSP1.2050*(SSP1.2050*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050= sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9, ely_total_SSP1_2020_q1= sum(cons_AC_SSP1.2020_q1*(SSP1.2020_q1*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050_q1= sum(cons_AC_SSP1.2050_q1*(SSP1.2050_q1*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050_q1= sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9, ely_total_SSP1_2020_q3= sum(cons_AC_SSP1.2020_q3*(SSP1.2020_q3*(pop_ssps_data_ssp1_2020)), na.rm = T)/1e9, ely_total_SSP1_2050_q3= sum(cons_AC_SSP1.2050_q3*(SSP1.2050_q3*(pop_ssps_data_ssp1_2050)), na.rm = T)/1e9, ely_total_SSP3_2050_q3= sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm = T)/1e9)

#

shape_ely_diff_gini <- dplyr::group_by(na.omit(shape_ely_diff_m), ISO3) %>% dplyr::reframe(ely_total_SSP2_2020= weighted.gini(cons_AC_SSP2.2020, (SSP2.2020*(pop_ssps_data_ssp2_2020))), ely_total_SSP2_2050= weighted.gini(cons_AC_SSP2.2050, (SSP2.2050*(pop_ssps_data_ssp2_2050))), ely_total_SSP5_2050= weighted.gini(cons_AC_SSP5.2050, (SSP5.2050*(pop_ssps_data_ssp5_2050))), ely_total_SSP1_2050= weighted.gini(cons_AC_SSP1.2050, (SSP1.2050*(pop_ssps_data_ssp1_2050))), ely_total_SSP3_2050= weighted.gini(cons_AC_SSP3.2050, (SSP3.2050*(pop_ssps_data_ssp3_2050))))
#

cols.num <- colnames(shape_ely_diff_gini)[2:6]

shape_ely_diff_gini[cols.num] <- sapply(shape_ely_diff_gini[cols.num],as.numeric)

shape_ely_diff_gini <- shape_ely_diff_gini %>% group_by(ISO3) %>% dplyr::summarise_at( .vars = colnames(.)[2:6], mean, na.rm=T)

colnames(shape_ely_diff_gini)[2:6] <- paste0(colnames(shape_ely_diff_gini)[2:6], "_gini")


data("wrld_simpl")
wrld_simpl <- subset(wrld_simpl, wrld_simpl$ISO3!="ATA")
wrld_simpl <- st_as_sf(wrld_simpl)

shape_ely_diff_s <- merge(shape_ely_diff_s, wrld_simpl, "ISO3")
shape_ely_diff_s <- st_as_sf(shape_ely_diff_s)

shape_ely_diff_s <- st_transform(shape_ely_diff_s, "ESRI:54009")
wrld_simpl <- st_transform(wrld_simpl, "ESRI:54009")
# 
# o2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ely_diff_s, aes(fill=ely_total_SSP2_2020), colour=NA, show.legend = F)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("AC-induced electricity consumption, 2020")
# 
# 
# a2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ely_diff_s, aes(fill=ely_total_SSP2_2050), colour=NA, show.legend = F)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("AC-induced electricity consumption, 2050, SSP2")
# 
# b2 <- ggplot()+
#   theme_void()+
#   geom_sf(data=shape_ely_diff_s, aes(fill=ely_total_SSP5_2050), colour=NA)+
#   geom_sf(data=wrld_simpl, colour="black", fill=NA, lwd=0.01)+
#   ggtitle("AC-induced electricity consumption, 2050, SSP5")
# 
# a_p2_2 <- o2 + a2 + b2 + plot_annotation(tag_levels = list(c('C', 'D', 'E'))) +  plot_layout(guides = "collect", ncol=1) &   scale_fill_binned(name="TWh/yr", type="viridis", trans="log10") & theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.width=unit(0.1,"npc"))
# 
# library(cowplot)
# 
# a_p2_m <- plot_grid(a_p2_1, a_p2_2)
# 
# ggsave("results/graphs_tables/maps_country.png", a_p2_m, scale=2)

#########

library(stargazer)

s_table <- shape_ely_diff_s %>% dplyr::select(ISO3, ely_total_SSP2_2020, ely_total_SSP2_2050, ely_total_SSP5_2050, ely_total_SSP2_2020_q1, ely_total_SSP2_2050_q1, ely_total_SSP5_2050_q1, ely_total_SSP2_2020_q3, ely_total_SSP2_2050_q3, ely_total_SSP5_2050_q3, ely_total_SSP1_2020, ely_total_SSP1_2050, ely_total_SSP3_2050, ely_total_SSP1_2020_q1, ely_total_SSP1_2050_q1, ely_total_SSP3_2050_q1, ely_total_SSP1_2020_q3, ely_total_SSP1_2050_q3, ely_total_SSP3_2050_q3) %>% mutate(ISO3=as.character(ISO3))
s_table$geometry <- NULL

s_table$ely_total_SSP2_2020 <- paste0(round(s_table$ely_total_SSP2_2020, 2), " (", round(s_table$ely_total_SSP2_2020_q1, 2), " - ", round(s_table$ely_total_SSP2_2020_q3, 2), ")")
s_table$ely_total_SSP2_2050 <- paste0(round(s_table$ely_total_SSP2_2050, 2), " (", round(s_table$ely_total_SSP2_2050_q1, 2), " - ", round(s_table$ely_total_SSP2_2050_q3, 2), ")")
s_table$ely_total_SSP5_2050 <- paste0(round(s_table$ely_total_SSP5_2050, 2), " (", round(s_table$ely_total_SSP5_2050_q1, 2), " - ", round(s_table$ely_total_SSP5_2050_q3, 2), ")")
s_table$ely_total_SSP1_2020 <- paste0(round(s_table$ely_total_SSP1_2020, 2), " (", round(s_table$ely_total_SSP1_2020_q1, 2), " - ", round(s_table$ely_total_SSP1_2020_q3, 2), ")")
s_table$ely_total_SSP1_2050 <- paste0(round(s_table$ely_total_SSP1_2050, 2), " (", round(s_table$ely_total_SSP1_2050_q1, 2), " - ", round(s_table$ely_total_SSP1_2050_q3, 2), ")")
s_table$ely_total_SSP3_2050 <- paste0(round(s_table$ely_total_SSP3_2050, 2), " (", round(s_table$ely_total_SSP3_2050_q1, 2), " - ", round(s_table$ely_total_SSP3_2050_q3, 2), ")")

s_table <- s_table %>% dplyr::select(-ely_total_SSP2_2020_q1, -ely_total_SSP2_2050_q1, -ely_total_SSP5_2050_q1, -ely_total_SSP2_2020_q3, -ely_total_SSP2_2050_q3, -ely_total_SSP5_2050_q3, -ely_total_SSP1_2020_q1, -ely_total_SSP1_2050_q1, -ely_total_SSP3_2050_q1, -ely_total_SSP1_2020_q3, -ely_total_SSP1_2050_q3, -ely_total_SSP3_2050_q3) %>% mutate(ISO3=as.character(ISO3))

###

s_table <- merge(s_table, shape_ely_diff_gini, by="ISO3")

s_table$ely_total_SSP1_2020 <- NULL
s_table$ely_total_SSP1_2020_gini <- NULL

s_table <- s_table[,c(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11)]

s_table <- s_table %>% dplyr::mutate_if(is.numeric, round, 2)

stargazer::stargazer(s_table, summary=F, out="results/graphs_tables/ELY_AC_consumpton_TWh.tex")

#############################

# shape_ac_s_top10 <- filter(shape_ac_s, shape_ac_s$POP2005 >=sort(shape_ac_s$POP2005, decreasing = T)[20])
# 
# shape_ac_s_top10$geometry <- NULL
# shape_ac_s_top10 <- dplyr::select(shape_ac_s_top10, SSP2.2020, SSP2.2050, SSP5.2050, ISO3)
# shape_ac_s_top10 <- reshape2::melt(shape_ac_s_top10, 4)
# 
# ggplot(shape_ac_s_top10)+
#   geom_col(aes(x=ISO3, y=value, fill=variable), position = "dodge")+
#   scale_y_continuous(labels=scales::label_percent())+
#   ylab("Projected AC adoption (%)")+
#   ggsci::scale_fill_npg(name="")
# 
# ggsave("results/graphs_tables/ac_col.png", width=10)
# 
# ##
# 
# shape_ely_diff_s_top10 <- filter(shape_ely_diff_s, shape_ely_diff_s$ely_total_SSP2_2050 >=sort(shape_ely_diff_s$ely_total_SSP2_2050, decreasing = T)[20])
# 
# shape_ely_diff_s_top10$geometry <- NULL
# shape_ely_diff_s_top10 <- dplyr::select(shape_ely_diff_s_top10, ely_total_SSP2_2010, ely_total_SSP2_2050, ely_total_SSP5_2050, ISO3)
# shape_ely_diff_s_top10 <- reshape2::melt(shape_ely_diff_s_top10, 4)
# 
# ggplot(shape_ely_diff_s_top10)+
#   geom_col(aes(x=ISO3, y=value, fill=gsub("ely_total_", "", variable)), position = "dodge")+
#   ggsci::scale_fill_npg(name="")+
#   ylab("AC-induced electr. cons. (TWh/yr.)")
# 
# ggsave("results/graphs_tables/ely_col.png", width=10)

##########

acglobtot = shape_ac %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2010_q1=weighted.mean(SSP2.2010_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q1=weighted.mean(SSP2.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q1=weighted.mean(SSP5.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2010_q3=weighted.mean(SSP2.2010_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2050_q3=weighted.mean(SSP2.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2050_q3=weighted.mean(SSP5.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T))

elyglobtot = dplyr::group_by(shape_ely_diff_m) %>% dplyr::summarise(ely_total_SSP2_2010= sum(cons_AC_SSP2.2010*(SSP2.2010*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050= sum(cons_AC_SSP2.2050*(SSP2.2050*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050= sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2010_q1= sum(cons_AC_SSP2.2010_q1*(SSP2.2010_q1*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q1= sum(cons_AC_SSP2.2050_q1*(SSP2.2050_q1*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q1= sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9, ely_total_SSP2_2010_q3= sum(cons_AC_SSP2.2010_q3*(SSP2.2010_q3*(pop_ssps_data_ssp2_2020)), na.rm = T)/1e9, ely_total_SSP2_2050_q3= sum(cons_AC_SSP2.2050_q3*(SSP2.2050_q3*(pop_ssps_data_ssp2_2050)), na.rm = T)/1e9, ely_total_SSP5_2050_q3= sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm = T)/1e9)

save("acglobtot", "elyglobtot", file=paste0(wd, "/results/glob_figures.Rdata"))

#### 

# pathways plot

paths <- dplyr::group_by(shape_ac, region) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2010_q1=weighted.mean(SSP2.2010_q1, pop_ssps_data_hist, na.rm=T), SSP2.2010_q3=weighted.mean(SSP2.2010_q3, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020_q1=weighted.mean(SSP2.2020_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020_q3=weighted.mean(SSP2.2020_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2030_q1=weighted.mean(SSP2.2030_q1, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2030_q3=weighted.mean(SSP2.2030_q3, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2040_q1=weighted.mean(SSP2.2040_q1, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2040_q3=weighted.mean(SSP2.2040_q3, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2050_q1=weighted.mean(SSP2.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2050_q3=weighted.mean(SSP2.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2010_q1=weighted.mean(SSP5.2010_q1, pop_ssps_data_hist, na.rm=T), SSP5.2010_q3=weighted.mean(SSP5.2010_q3, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2020_q1=weighted.mean(SSP5.2020_q1, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2020_q3=weighted.mean(SSP5.2020_q3, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2030_q1=weighted.mean(SSP5.2030_q1, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2030_q3=weighted.mean(SSP5.2030_q3, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2040_q1=weighted.mean(SSP5.2040_q1, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2040_q3=weighted.mean(SSP5.2040_q3, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T), SSP5.2050_q1=weighted.mean(SSP5.2050_q1, pop_ssps_data_ssp5_2050, na.rm=T), SSP5.2050_q3=weighted.mean(SSP5.2050_q3, pop_ssps_data_ssp5_2050, na.rm=T))


paths2 <- dplyr::group_by(shape_ac, region) %>% dplyr::summarise(SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2010_q1=weighted.mean(SSP1.2010_q1, pop_ssps_data_hist, na.rm=T), SSP1.2010_q3=weighted.mean(SSP1.2010_q3, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2020_q1=weighted.mean(SSP1.2020_q1, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2020_q3=weighted.mean(SSP1.2020_q3, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2030_q1=weighted.mean(SSP1.2030_q1, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2030_q3=weighted.mean(SSP1.2030_q3, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2040_q1=weighted.mean(SSP1.2040_q1, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2040_q3=weighted.mean(SSP1.2040_q3, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP1.2050_q1=weighted.mean(SSP1.2050_q1, pop_ssps_data_ssp1_2050, na.rm=T), SSP1.2050_q3=weighted.mean(SSP1.2050_q3, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2010_q1=weighted.mean(SSP3.2010_q1, pop_ssps_data_hist, na.rm=T), SSP3.2010_q3=weighted.mean(SSP3.2010_q3, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2020_q1=weighted.mean(SSP3.2020_q1, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2020_q3=weighted.mean(SSP3.2020_q3, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2030_q1=weighted.mean(SSP3.2030_q1, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2030_q3=weighted.mean(SSP3.2030_q3, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2040_q1=weighted.mean(SSP3.2040_q1, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2040_q3=weighted.mean(SSP3.2040_q3, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T), SSP3.2050_q1=weighted.mean(SSP3.2050_q1, pop_ssps_data_ssp3_2050, na.rm=T), SSP3.2050_q3=weighted.mean(SSP3.2050_q3, pop_ssps_data_ssp3_2050, na.rm=T))

paths <- merge(paths, paths2, by="region")

paths_glob <- dplyr::group_by(shape_ac) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2010_q1=weighted.mean(SSP2.2010_q1, pop_ssps_data_hist, na.rm=T), SSP2.2010_q3=weighted.mean(SSP2.2010_q3, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020_q1=weighted.mean(SSP2.2020_q1, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2020_q3=weighted.mean(SSP2.2020_q3, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2030_q1=weighted.mean(SSP2.2030_q1, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2030_q3=weighted.mean(SSP2.2030_q3, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2040_q1=weighted.mean(SSP2.2040_q1, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2040_q3=weighted.mean(SSP2.2040_q3, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2050_q1=weighted.mean(SSP2.2050_q1, pop_ssps_data_ssp2_2050, na.rm=T), SSP2.2050_q3=weighted.mean(SSP2.2050_q3, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2010_q1=weighted.mean(SSP5.2010_q1, pop_ssps_data_hist, na.rm=T), SSP5.2010_q3=weighted.mean(SSP5.2010_q3, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2020_q1=weighted.mean(SSP5.2020_q1, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2020_q3=weighted.mean(SSP5.2020_q3, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2030_q1=weighted.mean(SSP5.2030_q1, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2030_q3=weighted.mean(SSP5.2030_q3, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2040_q1=weighted.mean(SSP5.2040_q1, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2040_q3=weighted.mean(SSP5.2040_q3, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T), SSP5.2050_q1=weighted.mean(SSP5.2050_q1, pop_ssps_data_ssp5_2050, na.rm=T), SSP5.2050_q3=weighted.mean(SSP5.2050_q3, pop_ssps_data_ssp5_2050, na.rm=T))


paths_glob2 <- dplyr::group_by(shape_ac) %>% dplyr::summarise(SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2010_q1=weighted.mean(SSP1.2010_q1, pop_ssps_data_hist, na.rm=T), SSP1.2010_q3=weighted.mean(SSP1.2010_q3, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2020_q1=weighted.mean(SSP1.2020_q1, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2020_q3=weighted.mean(SSP1.2020_q3, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2030_q1=weighted.mean(SSP1.2030_q1, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2030_q3=weighted.mean(SSP1.2030_q3, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2040_q1=weighted.mean(SSP1.2040_q1, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2040_q3=weighted.mean(SSP1.2040_q3, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP1.2050_q1=weighted.mean(SSP1.2050_q1, pop_ssps_data_ssp1_2050, na.rm=T), SSP1.2050_q3=weighted.mean(SSP1.2050_q3, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2010_q1=weighted.mean(SSP3.2010_q1, pop_ssps_data_hist, na.rm=T), SSP3.2010_q3=weighted.mean(SSP3.2010_q3, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2020_q1=weighted.mean(SSP3.2020_q1, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2020_q3=weighted.mean(SSP3.2020_q3, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2030_q1=weighted.mean(SSP3.2030_q1, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2030_q3=weighted.mean(SSP3.2030_q3, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2040_q1=weighted.mean(SSP3.2040_q1, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2040_q3=weighted.mean(SSP3.2040_q3, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T), SSP3.2050_q1=weighted.mean(SSP3.2050_q1, pop_ssps_data_ssp3_2050, na.rm=T), SSP3.2050_q3=weighted.mean(SSP3.2050_q3, pop_ssps_data_ssp3_2050, na.rm=T))

paths_glob <- bind_cols(paths_glob, paths_glob2)

paths$type ="Region"
paths_glob$type ="Global"

paths <- bind_rows(paths, paths_glob)

paths <- reshape2::melt(paths, c(1, 62))
paths$ssp <- substr(paths$variable, 1, 4)
paths$ssp <- ifelse(paths$ssp=="SSP2", "SSP245", ifelse(paths$ssp=="SSP5","SSP585",  ifelse(paths$ssp=="SSP1", "SSP126", "SSP370")))

paths$year <- as.numeric(substr(paths$variable, 6, 9))

paths$q <- substr(paths$variable, 11, 12)
paths$q <- ifelse(paths$q=="", "q5", paths$q)

paths$ISO3 <- as.character(paths$region)
paths$ISO3 <- ifelse(is.na(paths$ISO3), " GLOBAL", paths$ISO3)

#####

paths_pivot <- pivot_wider(paths %>% dplyr::select(-variable), names_from = q, values_from = value)

#####

lines_a <- ggplot(paths_pivot)+
  theme_classic()+
  ggtitle("Projected evolution of residential AC penetration")+
  geom_line(aes(x=year, y=q5, group=ssp, colour=ssp, size=as.factor(desc(type))))+
  geom_ribbon(aes(x=year, ymin=q1, ymax=q3, group=ssp, fill=ssp), alpha=0.1)+
  facet_wrap(vars(ISO3), scales = "free", nrow=2)+
  scale_size_manual(values = c(.75, 1.5),guide = 'none')+
  scale_y_continuous(labels=scales::label_percent())+
  xlab("Year")+
  ylab("AC penetration rate")+
  scale_x_continuous(labels = c(2010, 2030, 2050), breaks = c(2010, 2030, 2050))+
  labs(caption = "Ribbon: IQR of CMIP6 GCMs") +
  theme(strip.background = element_blank()) 

#

shape_ely_diff <- merge(shape_ely_diff %>% dplyr::select(id, starts_with("cons_AC_SSP")), shape_ac %>% dplyr::select(id, region, ISO3, starts_with("SSP"), starts_with("pop_")), "id")

paths <- dplyr::group_by(shape_ely_diff, region) %>% dplyr::summarise(cons_AC_SSP2.2010=sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=sum(cons_AC_SSP2.2030*SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=sum(cons_AC_SSP2.2040*SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=sum(cons_AC_SSP5.2010*SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=sum(cons_AC_SSP5.2020*SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=sum(cons_AC_SSP5.2030*SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=sum(cons_AC_SSP5.2040*SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q1=sum(cons_AC_SSP2.2010_q1*SSP2.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q1=sum(cons_AC_SSP2.2020_q1*SSP2.2020_q1*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q1=sum(cons_AC_SSP2.2030_q1*SSP2.2030_q1*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q1=sum(cons_AC_SSP2.2040_q1*SSP2.2040_q1*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q1=sum(cons_AC_SSP2.2050_q1*SSP2.2050_q1*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q1=sum(cons_AC_SSP5.2010_q1*SSP5.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q1=sum(cons_AC_SSP5.2020_q1*SSP5.2020_q1*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q1=sum(cons_AC_SSP5.2030_q1*SSP5.2030_q1*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q1=sum(cons_AC_SSP5.2040_q1*SSP5.2040_q1*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q1=sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q3=sum(cons_AC_SSP2.2010_q3*SSP2.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q3=sum(cons_AC_SSP2.2020_q3*SSP2.2020_q3*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q3=sum(cons_AC_SSP2.2030_q3*SSP2.2030_q3*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q3=sum(cons_AC_SSP2.2040_q3*SSP2.2040_q3*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q3=sum(cons_AC_SSP2.2050_q3*SSP2.2050_q3*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q3=sum(cons_AC_SSP5.2010_q3*SSP5.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q3=sum(cons_AC_SSP5.2020_q3*SSP5.2020_q3*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q3=sum(cons_AC_SSP5.2030_q3*SSP5.2030_q3*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q3=sum(cons_AC_SSP5.2040_q3*SSP5.2040_q3*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q3=sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm=T))

paths2 <- dplyr::group_by(shape_ely_diff, region) %>% dplyr::summarise(cons_AC_SSP1.2010=sum(cons_AC_SSP1.2010*SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=sum(cons_AC_SSP1.2020*SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=sum(cons_AC_SSP1.2030*SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=sum(cons_AC_SSP1.2040*SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=sum(cons_AC_SSP3.2010*SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=sum(cons_AC_SSP3.2020*SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=sum(cons_AC_SSP3.2030*SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=sum(cons_AC_SSP3.2040*SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q1=sum(cons_AC_SSP1.2010_q1*SSP1.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q1=sum(cons_AC_SSP1.2020_q1*SSP1.2020_q1*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q1=sum(cons_AC_SSP1.2030_q1*SSP1.2030_q1*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q1=sum(cons_AC_SSP1.2040_q1*SSP1.2040_q1*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q1=sum(cons_AC_SSP1.2050_q1*SSP1.2050_q1*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q1=sum(cons_AC_SSP3.2010_q1*SSP3.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q1=sum(cons_AC_SSP3.2020_q1*SSP3.2020_q1*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q1=sum(cons_AC_SSP3.2030_q1*SSP3.2030_q1*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q1=sum(cons_AC_SSP3.2040_q1*SSP3.2040_q1*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q1=sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q3=sum(cons_AC_SSP1.2010_q3*SSP1.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q3=sum(cons_AC_SSP1.2020_q3*SSP1.2020_q3*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q3=sum(cons_AC_SSP1.2030_q3*SSP1.2030_q3*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q3=sum(cons_AC_SSP1.2040_q3*SSP1.2040_q3*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q3=sum(cons_AC_SSP1.2050_q3*SSP1.2050_q3*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q3=sum(cons_AC_SSP3.2010_q3*SSP3.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q3=sum(cons_AC_SSP3.2020_q3*SSP3.2020_q3*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q3=sum(cons_AC_SSP3.2030_q3*SSP3.2030_q3*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q3=sum(cons_AC_SSP3.2040_q3*SSP3.2040_q3*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q3=sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm=T))

paths <- merge(paths, paths2, by="region")

paths_glob <-dplyr::group_by(shape_ely_diff) %>% dplyr::summarise(cons_AC_SSP2.2010=sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=sum(cons_AC_SSP2.2030*SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=sum(cons_AC_SSP2.2040*SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=sum(cons_AC_SSP5.2010*SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=sum(cons_AC_SSP5.2020*SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=sum(cons_AC_SSP5.2030*SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=sum(cons_AC_SSP5.2040*SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q1=sum(cons_AC_SSP2.2010_q1*SSP2.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q1=sum(cons_AC_SSP2.2020_q1*SSP2.2020_q1*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q1=sum(cons_AC_SSP2.2030_q1*SSP2.2030_q1*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q1=sum(cons_AC_SSP2.2040_q1*SSP2.2040_q1*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q1=sum(cons_AC_SSP2.2050_q1*SSP2.2050_q1*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q1=sum(cons_AC_SSP5.2010_q1*SSP5.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q1=sum(cons_AC_SSP5.2020_q1*SSP5.2020_q1*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q1=sum(cons_AC_SSP5.2030_q1*SSP5.2030_q1*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q1=sum(cons_AC_SSP5.2040_q1*SSP5.2040_q1*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q1=sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q3=sum(cons_AC_SSP2.2010_q3*SSP2.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q3=sum(cons_AC_SSP2.2020_q3*SSP2.2020_q3*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q3=sum(cons_AC_SSP2.2030_q3*SSP2.2030_q3*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q3=sum(cons_AC_SSP2.2040_q3*SSP2.2040_q3*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q3=sum(cons_AC_SSP2.2050_q3*SSP2.2050_q3*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q3=sum(cons_AC_SSP5.2010_q3*SSP5.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q3=sum(cons_AC_SSP5.2020_q3*SSP5.2020_q3*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q3=sum(cons_AC_SSP5.2030_q3*SSP5.2030_q3*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q3=sum(cons_AC_SSP5.2040_q3*SSP5.2040_q3*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q3=sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm=T)) 

paths_glob2 <-dplyr::group_by(shape_ely_diff) %>% dplyr::summarise(cons_AC_SSP1.2010=sum(cons_AC_SSP1.2010*SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=sum(cons_AC_SSP1.2020*SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=sum(cons_AC_SSP1.2030*SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=sum(cons_AC_SSP1.2040*SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=sum(cons_AC_SSP3.2010*SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=sum(cons_AC_SSP3.2020*SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=sum(cons_AC_SSP3.2030*SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=sum(cons_AC_SSP3.2040*SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q1=sum(cons_AC_SSP1.2010_q1*SSP1.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q1=sum(cons_AC_SSP1.2020_q1*SSP1.2020_q1*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q1=sum(cons_AC_SSP1.2030_q1*SSP1.2030_q1*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q1=sum(cons_AC_SSP1.2040_q1*SSP1.2040_q1*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q1=sum(cons_AC_SSP1.2050_q1*SSP1.2050_q1*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q1=sum(cons_AC_SSP3.2010_q1*SSP3.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q1=sum(cons_AC_SSP3.2020_q1*SSP3.2020_q1*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q1=sum(cons_AC_SSP3.2030_q1*SSP3.2030_q1*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q1=sum(cons_AC_SSP3.2040_q1*SSP3.2040_q1*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q1=sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q3=sum(cons_AC_SSP1.2010_q3*SSP1.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q3=sum(cons_AC_SSP1.2020_q3*SSP1.2020_q3*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q3=sum(cons_AC_SSP1.2030_q3*SSP1.2030_q3*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q3=sum(cons_AC_SSP1.2040_q3*SSP1.2040_q3*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q3=sum(cons_AC_SSP1.2050_q3*SSP1.2050_q3*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q3=sum(cons_AC_SSP3.2010_q3*SSP3.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q3=sum(cons_AC_SSP3.2020_q3*SSP3.2020_q3*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q3=sum(cons_AC_SSP3.2030_q3*SSP3.2030_q3*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q3=sum(cons_AC_SSP3.2040_q3*SSP3.2040_q3*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q3=sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm=T)) 

paths_glob <- bind_cols(paths_glob, paths_glob2)

paths$type ="Region"
paths_glob$type ="Global"

paths <- bind_rows(paths, paths_glob)

paths <- reshape2::melt(paths, c(1, 62))
paths$variable <- gsub("cons_AC_", "", paths$variable)
paths$ssp <- substr(paths$variable, 1, 4)
paths$ssp <- ifelse(paths$ssp=="SSP2", "SSP245", ifelse(paths$ssp=="SSP5","SSP585",  ifelse(paths$ssp=="SSP1", "SSP126", "SSP370")))
paths$year <- as.numeric(substr(paths$variable, 6, 9))

paths$q <- substr(paths$variable, 11, 12)
paths$q <- ifelse(paths$q=="", "q5", paths$q)

paths$ISO3 <- as.character(paths$region)
paths$ISO3 <- ifelse(is.na(paths$ISO3), " GLOBAL", paths$ISO3)

paths_pivot <- pivot_wider(paths %>% dplyr::select(-variable), names_from = q, values_from = value)

writexl::write_xlsx(paths_pivot, "paths_pivot.xlsx")

lines_b <- ggplot(paths_pivot)+
  theme_classic()+
  ggtitle("Projected evolution of residential AC electricity consumption")+
  geom_line(aes(x=year, y=q5/1e9, group=ssp, colour=ssp, size=as.factor(desc(type))), show.legend = F)+
  geom_ribbon(aes(x=year, ymin=q1/1e9, ymax=q3/1e9, group=ssp, fill=ssp), alpha=0.1)+
  facet_wrap(vars(ISO3), scales = "free", nrow=2)+
  scale_size_manual(values = c(0.75, 1.5),guide = 'none')+
  xlab("Year")+
  ylab("AC electricity consumption (TWh/yr.)")+
  scale_x_continuous(labels = c(2010, 2030, 2050), breaks = c(2010, 2030, 2050))+
  labs(caption = "Ribbon: IQR of CMIP6 GCMs") +
  theme(strip.background = element_blank()) 

library(patchwork)

(lines_a + scale_colour_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404")) + scale_fill_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))) + (lines_b + scale_colour_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404")) + scale_fill_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404"))) + plot_layout(guides = "collect", ncol = 1) + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", legend.direction = "horizontal")  & guides(colour = guide_legend(nrow = 1))

ggsave("results/graphs_tables/lines_plot.pdf", scale=2, height = 5, width = 4.5)

#

emissions_figure <- read_rds("results/emissions_figure.rds")


(lines_b + scale_colour_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404")) + scale_fill_manual(name="Scenario", values=c("#fcfc65", "#facf96", "#e38202", "#7d0404")) + theme(legend.position = "none")) + (emissions_figure + ggtitle("Projected evolution of GHG emissions from residential AC") + ylab("AC GHG emissions (Mt CO2e/yr.)") + theme(legend.position = "bottom", legend.direction = "horizontal") + guides(colour = guide_legend(nrow = 1))) + plot_layout(ncol=1) + plot_annotation(tag_levels = "A")

ggsave("results/graphs_tables/lines_plot_new.pdf", scale=2, height = 5, width = 4.5)


#

paths_si <- dplyr::group_by(shape_ely_diff, ISO3) %>% dplyr::summarise(cons_AC_SSP2.2010=sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=sum(cons_AC_SSP2.2030*SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=sum(cons_AC_SSP2.2040*SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=sum(cons_AC_SSP5.2010*SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=sum(cons_AC_SSP5.2020*SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=sum(cons_AC_SSP5.2030*SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=sum(cons_AC_SSP5.2040*SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T))

paths_si2 <- dplyr::group_by(shape_ely_diff, ISO3) %>% dplyr::summarise(cons_AC_SSP1.2010=sum(cons_AC_SSP1.2010*SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=sum(cons_AC_SSP1.2020*SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=sum(cons_AC_SSP1.2030*SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=sum(cons_AC_SSP1.2040*SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=sum(cons_AC_SSP3.2010*SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=sum(cons_AC_SSP3.2020*SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=sum(cons_AC_SSP3.2030*SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=sum(cons_AC_SSP3.2040*SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T))

paths_si <- merge(paths_si, paths_si2, by="ISO3")

paths_si <- reshape2::melt(paths_si, c(1))
paths_si$variable <- gsub("cons_AC_", "", paths_si$variable)
paths_si$ssp <- substr(paths_si$variable, 1, 4)
paths_si$year <- as.numeric(substr(paths_si$variable, 6, 9))
paths_si$variable <- NULL

paths_si_ely <- paths_si

###

paths_si <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T))

paths_si2 <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

paths_si <- merge(paths_si, paths_si2, by="ISO3")

paths_si <- reshape2::melt(paths_si, c(1))
paths_si$ssp <- substr(paths_si$variable, 1, 4)
paths_si$year <- as.numeric(substr(paths_si$variable, 6, 9))
paths_si$variable <- NULL

paths_si_ac <- paths_si

paths_si <- paths_si_ely
paths_si$value_ac <- paths_si_ac$value

###

library(modelsummary)

paths_si$value <- paths_si$value/1e9
paths_si$value_ac <- paths_si$value_ac*100

paths_si <- paths_si %>% 
  mutate_if(is.numeric, round, 1)

emptycol <- function(x) " "

paths_si <- reshape2::melt(paths_si, c(1,3,4))

paths_si$ISO3 <- as.character(paths_si$ISO3)

paths_si <- na.omit(paths_si)

paths_si$variable <- as.character(paths_si$variable)
paths_si$variable[paths_si$variable=="value"] <- "Electricity consumption for AC (TWh/yr.)"
paths_si$variable[paths_si$variable=="value_ac"] <- "AC penetration (% of HHs)"

paths_si <- filter(paths_si, year==2050)

datasummary( value * ISO3 ~ variable * ssp * (Mean), data = paths_si, longtable = TRUE, output = "results/graphs_tables/si_ctry_table.tex")

#############

#faceted by region, ac penetration by decile of income (today and in two scenarios)

shape_ac <- shape_ac %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

paths <- dplyr::group_by(shape_ac, region, decile) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T))

paths_glob <- dplyr::group_by(shape_ac, decile) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T), SSP2.2020=weighted.mean(SSP2.2020, pop_ssps_data_ssp2_2020, na.rm=T), SSP2.2030=weighted.mean(SSP2.2030, pop_ssps_data_ssp2_2030, na.rm=T), SSP2.2040=weighted.mean(SSP2.2040, pop_ssps_data_ssp2_2040, na.rm=T), SSP2.2050=weighted.mean(SSP2.2050, pop_ssps_data_ssp2_2050, na.rm=T), SSP5.2010=weighted.mean(SSP5.2010, pop_ssps_data_hist, na.rm=T), SSP5.2020=weighted.mean(SSP5.2020, pop_ssps_data_ssp5_2020, na.rm=T), SSP5.2030=weighted.mean(SSP5.2030, pop_ssps_data_ssp5_2030, na.rm=T), SSP5.2040=weighted.mean(SSP5.2040, pop_ssps_data_ssp5_2040, na.rm=T), SSP5.2050=weighted.mean(SSP5.2050, pop_ssps_data_ssp5_2050, na.rm=T))

###

paths2 <- dplyr::group_by(shape_ac, region, decile) %>% dplyr::summarise(SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

paths_glob2 <- dplyr::group_by(shape_ac, decile) %>% dplyr::summarise(SSP1.2010=weighted.mean(SSP1.2010, pop_ssps_data_hist, na.rm=T), SSP1.2020=weighted.mean(SSP1.2020, pop_ssps_data_ssp1_2020, na.rm=T), SSP1.2030=weighted.mean(SSP1.2030, pop_ssps_data_ssp1_2030, na.rm=T), SSP1.2040=weighted.mean(SSP1.2040, pop_ssps_data_ssp1_2040, na.rm=T), SSP1.2050=weighted.mean(SSP1.2050, pop_ssps_data_ssp1_2050, na.rm=T), SSP3.2010=weighted.mean(SSP3.2010, pop_ssps_data_hist, na.rm=T), SSP3.2020=weighted.mean(SSP3.2020, pop_ssps_data_ssp3_2020, na.rm=T), SSP3.2030=weighted.mean(SSP3.2030, pop_ssps_data_ssp3_2030, na.rm=T), SSP3.2040=weighted.mean(SSP3.2040, pop_ssps_data_ssp3_2040, na.rm=T), SSP3.2050=weighted.mean(SSP3.2050, pop_ssps_data_ssp3_2050, na.rm=T))

paths <- merge(paths, paths2, by=c("region", "decile"))
paths_glob <- merge(paths_glob, paths_glob2, by="decile")


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

lines_a <- ggplot(paths)+
  theme_classic()+
  geom_col(aes(x=level, y=(value)*100, group=decile, fill=decile), size=1.5, colour="lightgrey", lwd=0.01,width=0.8,    
           position=position_dodge(.8))+
  facet_wrap(vars(ISO3), ncol=4)+
  scale_fill_brewer(name="Income quintile", palette="Blues")+
  xlab("Scenario")+
  ylab("AC penetration rate")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region.")


lines_a <- ggplot(paths)+
  theme_classic()+
  geom_line(aes(x=decile, y=(value)*100, group=level, colour=level), size=.5)+
  geom_point(aes(x=decile, y=(value)*100, group=level, colour=level), size=1.5)+
  facet_wrap(vars(ISO3), ncol=4)+
  scale_colour_manual(name="", values=c("grey", "#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  xlab("Income quintile")+
  ylab("AC penetration rate (%)")+
  theme(legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region.")

write_rds(paths, "results/quntiles.rds")

###

#faceted by region, ely penetration by decile of income (today and in two scenarios)

shape_ely_diff_m$region <- shape_ely_diff_m$region.x

shape_ely_diff <- shape_ely_diff_m %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

#shape_ely_diff <- shape_ely_diff_m %>% group_by(region) %>% mutate(decile=ntile(`CDD_126_2015_CMCC-CM2-SR5`, 5))

shape_ely_diff <- na.omit(shape_ely_diff)

shape_ely_diff$decile <- as.factor(shape_ely_diff$decile)

paths <- dplyr::group_by(shape_ely_diff, region, decile) %>% dplyr::summarise(cons_AC_SSP2.2010=weighted.mean(cons_AC_SSP2.2010, SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=weighted.mean(cons_AC_SSP2.2020, SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=weighted.mean(cons_AC_SSP2.2030, SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=weighted.mean(cons_AC_SSP2.2040, SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=weighted.mean(cons_AC_SSP2.2050, SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=weighted.mean(cons_AC_SSP5.2010, SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=weighted.mean(cons_AC_SSP5.2020, SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=weighted.mean(cons_AC_SSP5.2030, SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=weighted.mean(cons_AC_SSP5.2040, SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=weighted.mean(cons_AC_SSP5.2050, SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP1.2010=weighted.mean(cons_AC_SSP1.2010, SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=weighted.mean(cons_AC_SSP1.2020, SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=weighted.mean(cons_AC_SSP1.2030, SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=weighted.mean(cons_AC_SSP1.2040, SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=weighted.mean(cons_AC_SSP1.2050, SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=weighted.mean(cons_AC_SSP3.2010, SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=weighted.mean(cons_AC_SSP3.2020, SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=weighted.mean(cons_AC_SSP3.2030, SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=weighted.mean(cons_AC_SSP3.2040, SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=weighted.mean(cons_AC_SSP3.2050, SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T))

paths_glob <- dplyr::group_by(shape_ely_diff, decile) %>% dplyr::summarise(cons_AC_SSP2.2010=weighted.mean(cons_AC_SSP2.2010, SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=weighted.mean(cons_AC_SSP2.2020, SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=weighted.mean(cons_AC_SSP2.2030, SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=weighted.mean(cons_AC_SSP2.2040, SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=weighted.mean(cons_AC_SSP2.2050, SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=weighted.mean(cons_AC_SSP5.2010, SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=weighted.mean(cons_AC_SSP5.2020, SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=weighted.mean(cons_AC_SSP5.2030, SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=weighted.mean(cons_AC_SSP5.2040, SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=weighted.mean(cons_AC_SSP5.2050, SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP1.2010=weighted.mean(cons_AC_SSP1.2010, SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=weighted.mean(cons_AC_SSP1.2020, SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=weighted.mean(cons_AC_SSP1.2030, SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=weighted.mean(cons_AC_SSP1.2040, SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=weighted.mean(cons_AC_SSP1.2050, SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=weighted.mean(cons_AC_SSP3.2010, SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=weighted.mean(cons_AC_SSP3.2020, SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=weighted.mean(cons_AC_SSP3.2030, SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=weighted.mean(cons_AC_SSP3.2040, SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=weighted.mean(cons_AC_SSP3.2050, SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T))


paths$type ="Region"
paths_glob$type ="Global"

paths <- bind_rows(paths, paths_glob)

paths <- reshape2::melt(paths, c(1, 2, 23))
paths$ssp <- substr(paths$variable, 9, 12)
paths$year <- as.numeric(substr(paths$variable, 14, 17))

paths$ISO3 <- as.character(paths$region)
paths$ISO3 <- ifelse(is.na(paths$ISO3), " GLOBAL", paths$ISO3)

paths$decile <- as.factor(paths$decile)

paths <- filter(paths, year==2020 | year==2050)

paths$level <- ifelse(paths$ssp=="SSP2" & paths$year==2020, "2020", ifelse(paths$ssp=="SSP2" & paths$year==2050, "2050, SSP245",  ifelse(paths$ssp=="SSP5" & paths$year==2050, "2050, SSP585", ifelse(paths$ssp=="SSP1" & paths$year==2050, "2050, SSP126", ifelse(paths$ssp=="SSP3" & paths$year==2050, "2050, SSP370", NA)))))

paths <- filter(paths, !is.na(level))
paths$level <- as.factor(paths$level)

lines_b <- ggplot(paths)+
  theme_classic()+
  geom_col(aes(x=level, y=(value), group=decile, fill=decile), size=1.5, colour="lightgrey", lwd=0.01,width=0.8,    
           position=position_dodge(.8))+
  facet_wrap(vars(ISO3), ncol=4)+
  scale_fill_brewer(name="Income quintile", palette="Blues")+
  xlab("Scenario")+
  ylab("Mean AC electricity consumption per HH")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region.")

lines_b <- ggplot(paths)+
  theme_classic()+
  geom_line(aes(x=decile, y=(value), group=level, colour=level), size=.5)+
  geom_point(aes(x=decile, y=(value), group=level, colour=level), size=1.5)+
  facet_wrap(vars(ISO3), ncol=4)+
  scale_colour_manual(name="", values=c("grey", "#fcfc65", "#facf96", "#e38202", "#7d0404"))+
  xlab("Income quintile")+
  ylab("Mean AC electricity consumption (kWh/hh./yr.)")+
  theme(legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region.")

write_rds(paths, "results/quntiles_ely.rds")

library(patchwork)

(lines_a + lines_b) + plot_layout(guides = "collect", ncol = 1) + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", legend.direction = "horizontal")  & guides(colour = guide_legend(nrow = 1))

# ggsave("results/graphs_tables/quntiles_plot.pdf", scale=1.25, height = 8, width = 6)


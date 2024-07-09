# table drivers

## 1) Load libraries and data ##
library(tidyverse)
library(reshape2)
library(sf)
#library(gg.layers)
####

setwd(wd)

load("supporting_data/data_for_global_spline_v2.Rds")

############

shape <- st_as_sf(shape)
shape <- filter(shape, !is.na(region))

shape_gdp <- dplyr::select(shape, contains("GDP"), region)
shape_gdp$geometry <- NULL

shape_gdp <- melt(shape_gdp, 26)
shape_gdp$ssp <- substr(shape_gdp$variable, 5, 8)
shape_gdp$year <- substr(shape_gdp$variable, 10, 13)

shape_gdp <- filter(shape_gdp, year=="2020" | year=="2050")

###############

shape_pop <- dplyr::select(shape, contains("pop"), region)
shape_pop$POP2005 = NULL
shape_pop$geometry <- NULL

shape_pop <- melt(shape_pop, 27)
shape_pop$ssp <- gsub("pop_ssps_data_", "", shape_pop$variable)
shape_pop$ssp <- gsub("_2020|_2030|_2040|_2050", "", shape_pop$ssp)
shape_pop <- filter(shape_pop, ssp!="hist")

shape_pop$year <- gsub("pop_ssps_data_ssp1_", "", shape_pop$variable)
shape_pop$year <- gsub("pop_ssps_data_ssp4_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp3_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp2_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp5_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_", "", shape_pop$year)

shape_pop <- filter(shape_pop, year=="2020" | year=="2050")

#

# calculate IQR of GCMs

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

shape <- ungroup(shape)
shape$geometry <- NULL

future_acc <- data.frame(region=shape$region)

for(x in c("126", "245", "370", "585")){
  for(y in seq(2015, 2050, 1)){
    
    col1_columns <- grepl(y, names(shape)) & grepl(x, names(shape)) & (names(shape)!=paste0(x, ".", y)) & grepl("CDD", names(shape))
    
    varname <- paste0("CDD_", x, "_", y)
    
    future_acc <-bind_cols(future_acc, shape %>% dplyr::select(colnames(shape)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.5)))
    
  }}


shape_cdd <- future_acc

shape_cdd <- melt(shape_cdd, 1)
shape_cdd$ssp <- substr(shape_cdd$variable, 5, 7)
shape_cdd$year <- substr(shape_cdd$variable, 9, 12)

shape_cdd <- filter(shape_cdd, year=="2020" | year=="2050")

shape_cdd$id <- NULL

####

shape <- ungroup(shape)

future_acc <- data.frame(region=shape$region)

for(x in c("126", "245", "370", "585")){
  for(y in seq(2015, 2050, 1)){
    
    col1_columns <- grepl(y, names(shape)) & grepl(x, names(shape)) & (names(shape)!=paste0(x, ".", y)) & grepl("HDD", names(shape))
    
    varname <- paste0("HDD_", x, "_", y)
    
    future_acc <-bind_cols(future_acc, shape %>% dplyr::select(colnames(shape)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.5)))
    
  }}


shape_hdd <- future_acc

shape_hdd <- melt(shape_hdd, 1)
shape_hdd$ssp <- substr(shape_hdd$variable, 5, 7)
shape_hdd$year <- substr(shape_hdd$variable, 9, 12)

shape_hdd <- filter(shape_hdd, year=="2020" | year=="2050")

shape_hdd$id <- NULL

shape_hdd <- filter(shape_hdd, !is.na(region))

#

shape_age <- dplyr::select(shape, contains("age"), region)
shape_age$geometry <- NULL

shape_age <- melt(shape_age, 26)
shape_age$ssp <- substr(shape_age$variable, 5, 8)
shape_age$year <- substr(shape_age$variable, 10, 13)

shape_age <- filter(shape_age, year=="2020" | year=="2050")

shape_age <- filter(shape_age, !is.na(region))

shape_age <- shape_age %>% 
  group_by(region, ssp, year) %>%
  mutate(med_age = median(as.numeric(value), na.rm = T))

#

shape_edu <- dplyr::select(shape, contains("edu"), region)
shape_edu$geometry <- NULL

shape_edu <- melt(shape_edu, 26)
shape_edu$ssp <- substr(shape_edu$variable, 5, 8)
shape_edu$year <- substr(shape_edu$variable, 10, 13)

shape_edu <- filter(shape_edu, year=="2020" | year=="2050")

shape_edu <- filter(shape_edu, !is.na(region))

shape_edu <- shape_edu %>% 
  group_by(region, ssp, year) %>%
  mutate(med_edu = median(as.numeric(value), na.rm = T))

#################
#################


shape_age = shape_age %>% group_by(region, ssp, year) %>% dplyr::summarise(value=mean(value, na.rm=T), variable="Age")

shape_edu = shape_edu %>% group_by(region, ssp, year) %>% dplyr::summarise(value=mean(value, na.rm=T), variable= "Education")

shape_cdd = shape_cdd %>% group_by(region, ssp, year) %>% dplyr::summarise(value=mean(value, na.rm=T), variable= "CDDs")

shape_hdd = shape_hdd %>% group_by(region, ssp, year) %>% dplyr::summarise(value=mean(value, na.rm=T), variable= "HDDs")

shape_gdp = shape_gdp %>% group_by(region, ssp, year) %>% dplyr::summarise(value=median(value, na.rm=T), variable= "Per-capita GDP")

shape_pop = shape_pop %>% group_by(region, ssp, year) %>% dplyr::summarise(value=sum(value, na.rm=T), variable= "Population")

####

shape_m = bind_rows(shape_age, shape_cdd, shape_edu, shape_gdp, shape_hdd, shape_pop)

shape_m$ssp_full = ifelse(shape_m$ssp=="ssp1", "SSP1", ifelse(shape_m$ssp=="ssp2", "SSP2", ifelse(shape_m$ssp=="ssp3", "SSP3", ifelse(shape_m$ssp=="ssp4", "SSP4", ifelse(shape_m$ssp=="ssp5", "SSP5", ifelse(shape_m$ssp=="126", "SSP1", ifelse(shape_m$ssp=="245", "SSP2", ifelse(shape_m$ssp=="370", "SSP3", ifelse(shape_m$ssp=="585", "SSP5", shape_m$ssp)))))))))

shape_m$ssp <- shape_m$ssp_full
shape_m$ssp_full <- NULL

library(modelsummary)
 
shape_m$value <- ifelse(shape_m$variable=="Population", shape_m$value/1e9, shape_m$value)
shape_m$value <- ifelse(shape_m$variable=="Per-capita GDP", exp(shape_m$value), shape_m$value)

shape_m <- shape_m %>% 
  mutate_if(is.numeric, round, 1)

shape_m <- filter(shape_m, ssp!="SSP4")

datasummary(value * region * year * ssp ~  variable * (Mean),
            data = shape_m, output="results/drivers_evolution.tex")

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
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

setwd(wd)

##############

# extensive margin

###

load(paste0(wd, "/results/xgboost_models_benchmarks_jan24.Rdata"))
load(paste0(wd, "/results/xgboost_models_jan24.Rdata"))
load(paste0(wd, "/supporting_data/data_for_global_spline_v2.Rds"))

###

# check NAs

# library(pbapply)
# 
# l <- pblapply(5:(ncol(shape)-1), function(L) tapply(as.data.frame(shape[,c(L)] %>% st_set_geometry(NULL))[,1], shape$ISO3, function(X){
#   a <- as.numeric((sum(is.na(X))))
#   }))
# l <- as.data.frame( do.call(cbind, l))
# l[l==0] <- NA
# colnames(l) <- colnames(shape)[5:(ncol(shape)-1)]
# l <- l[!complete.cases(l)]

###

adj_gtap <- read.csv("adj_gtap/exp_gdp_capita_markups.csv")

adj_gtap$REG <- toupper(adj_gtap$REG)

shape <- merge(shape, adj_gtap, by.x="ISO3", by.y="REG", all.x=T)

shape <- shape %>% group_by(region) %>% mutate(markup = ifelse(is.na(markup), median(markup, na.rm=T), markup))

#shape <- mutate_at(shape, vars(starts_with("GDP_")), ~log(exp(.)*markup))

shape$markup <- NULL

shape$macroregion <- countrycode::countrycode(shape$ISO3, 'iso3c', 'continent')

# shape <- mutate_at(shape, vars(contains('GDP')), funs(ifelse(exp(.)>0.75e5, log(0.75e5), .)))
# shape <- mutate_at(shape, vars(contains('CDD')), funs(log((.*100)+1)))
# shape <- mutate_at(shape, vars(contains('CDD')), funs(ifelse(.>7.309, 7.309, .)))
# 

library(imputeTS)
library(zoo)

shape$X <- NULL

shape <- shape %>% 
  group_by(ISO3) %>% 
  mutate_if(is.numeric, zoo::na.aggregate)

shape <- shape %>% 
  group_by(region) %>% 
  mutate_if(is.numeric, zoo::na.aggregate)

shape <- na.omit(shape)

shape_all <- shape

shape$geometry <- NULL

###

# pops <- read_rds("pop_ssps_data_hist.Rds")
# 
# shape <- merge(shape, pops, "id")

#

orig_data_bk <- shape

##########

hhsize = read.csv("supporting_data/hhsize.csv", sep = ";", encoding = "Latin-1", stringsAsFactors = F)
hhsize$year <- as.Date(hhsize$year, format = "%d/%M/%Y")
hhsize$year <- lubridate::year(hhsize$year)
hhsize = hhsize %>% group_by(country) %>% filter(year==max(year)) %>% ungroup()
hhsize$value <- gsub(",", ".", hhsize$value)
hhsize$value <- as.numeric(hhsize$value)
hhsize$Country = countrycode::countrycode(hhsize$country, 'country.name', 'iso3c')
hhsize <- dplyr::select(hhsize, Country, value)
hhsize <- as.data.frame(hhsize)

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)

wrld_simpl_hhsize = base::merge(wrld_simpl, hhsize, by.x="ISO3", by.y="Country", all.x=T)
wrld_simpl_hhsize$value = ifelse(is.na(wrld_simpl_hhsize$value), mean(wrld_simpl_hhsize$value, na.rm=T), wrld_simpl_hhsize$value)

pop_ssps <- list.files(path="supporting_data/pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- lapply(pop_ssps, stack)

wrld_simpl_hhsize = fasterize::fasterize(st_as_sf(wrld_simpl_hhsize), pop_ssps_data[[2]][[1]], field="value")

shape$hhsize = exactextractr::exact_extract(wrld_simpl_hhsize, shape_all$geometry, "mean")

shape <- as.data.frame(shape)

###########

output3 <- list()

# loop for all ssps and time-steps

load(paste0(wd, "/supporting_data/models_list.Rds"))
cmip6_models <- models_list

####

# source("decomposition_drivers.R")

####

for (model in cmip6_models){
print(model)
  
  output <- list()
  
for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  orig_data_bk <- shape
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    #orig_data$weight = 0
    
    orig_data$ln_total_exp_usd_2011 = shape[,paste0("GDP_", ssp, "_", year)]
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, year), "_", model)]
    
    orig_data$ely_p_usd_2011 <- shape$elyprc
    
    orig_data$inc_q <- cut(exp(orig_data$ln_total_exp_usd_2011), terc_inc, include.lowest = T, labels=c(1,2,3))
    
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)<terc_inc[1], 1, orig_data$inc_q)
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)>terc_inc[4], 3, orig_data$inc_q)
    
    orig_data$inc_q <- as.factor(orig_data$inc_q)
    
    orig_data$cdd_q <- cut(exp(orig_data$mean_CDD18_db), terc_cdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)<terc_cdd[1], 1, orig_data$cdd_q)
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)>terc_cdd[4], 3, orig_data$cdd_q)
    
    orig_data$cdd_q <- as.factor(orig_data$cdd_q)
    
    orig_data$hdd_q <-  cut(exp(orig_data$mean_HDD18_db), terc_hdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)<terc_hdd[1], 1, orig_data$hdd_q)
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)>terc_hdd[4], 3, orig_data$hdd_q)
    
    orig_data$hdd_q <- as.factor(orig_data$hdd_q)
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", year)] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", (year))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}
  output3[[as.character(model)]] <- output
  
  }

##################
##################

future_ac_adoption <- do.call(cbind, lapply(output3, as.data.frame))

future_ac_adoption[future_ac_adoption<=0] <- NA

shape_ac <- bind_cols(shape, future_ac_adoption)

##############

# intensive margin

shape_ac <- na_mean(shape_ac)

orig_data <- shape_ac
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

for (model in cmip6_models){
  print(model)
  output <- list()

# loop for all ssps and time-steps

  for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
    
    print(ssp)
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = shape_ac[,paste0(model, ".", ssp, ".", year)] 
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", year)])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, year), "_", model)]
    
    orig_data$ely_p_usd_2011 <- shape$elyprc
    
    orig_data$inc_q <- cut(exp(orig_data$ln_total_exp_usd_2011), terc_inc, include.lowest = T, labels=c(1,2,3))
    
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)<terc_inc[1], 1, orig_data$inc_q)
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)>terc_inc[4], 3, orig_data$inc_q)
    
    orig_data$inc_q <- as.factor(orig_data$inc_q)
    
    orig_data$cdd_q <- cut(exp(orig_data$mean_CDD18_db), terc_cdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)<terc_cdd[1], 1, orig_data$cdd_q)
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)>terc_cdd[4], 3, orig_data$cdd_q)
    
    orig_data$cdd_q <- as.factor(orig_data$cdd_q)
    
    orig_data$hdd_q <-  cut(exp(orig_data$mean_HDD18_db), terc_hdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)<terc_hdd[1], 1, orig_data$hdd_q)
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)>terc_hdd[4], 3, orig_data$hdd_q)
    
    orig_data$hdd_q <- as.factor(orig_data$hdd_q)
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", year)] 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", year)] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", (year))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}
  
  output3[[as.character(model)]] <- output
  
}

future_ac_cons <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons <- exp(future_ac_cons)
future_ac_cons[future_ac_cons<=0] <- NA

# future_ac_cons <- future_ac_cons %>%
#   mutate_if(is.numeric, ~ replace(., .>10000, 10000))

#write_rds(future_ac_cons, "future_ac_cons.rds")

shape_ely <- bind_cols(shape_all %>% st_set_geometry(NULL), future_ac_cons)

r <- raster(ncol=4320, nrow=2160); r[] <- 1:ncell(r)

shape_prova <- fasterize(shape_ely %>% st_set_geometry(shape_all$geometry), r, "CMCC-CM2-SR5.SSP2.2020")

data("wrld_simpl")
plot(shape_prova)
plot(wrld_simpl, add=T)

#

orig_data <- shape_ac
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

####

output3 <- list()

# loop for all ssps and time-steps

for (model in cmip6_models){
  print(model)
  output <- list()

# loop for all ssps and time-steps

  for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
    
    print(ssp)
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
    
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 0
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", year)])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, year), "_", model)] +1)

    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, year), "_", model)]
    
    orig_data$ely_p_usd_2011 <- shape$elyprc
    
    orig_data$inc_q <- cut(exp(orig_data$ln_total_exp_usd_2011), terc_inc, include.lowest = T, labels=c(1,2,3))
    
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)<terc_inc[1], 1, orig_data$inc_q)
    orig_data$inc_q <- ifelse(exp(orig_data$ln_total_exp_usd_2011)>terc_inc[4], 3, orig_data$inc_q)
    
    orig_data$inc_q <- as.factor(orig_data$inc_q)
    
    orig_data$cdd_q <- cut(exp(orig_data$mean_CDD18_db), terc_cdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)<terc_cdd[1], 1, orig_data$cdd_q)
    orig_data$cdd_q <- ifelse(exp(orig_data$mean_CDD18_db)>terc_cdd[4], 3, orig_data$cdd_q)
    
    orig_data$cdd_q <- as.factor(orig_data$cdd_q)
    
    orig_data$hdd_q <-  cut(exp(orig_data$mean_HDD18_db), terc_hdd, include.lowest = T, labels=c(1,2,3))
    
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)<terc_hdd[1], 1, orig_data$hdd_q)
    orig_data$hdd_q <- ifelse(exp(orig_data$mean_HDD18_db)>terc_hdd[4], 3, orig_data$hdd_q)
    
    orig_data$hdd_q <- as.factor(orig_data$hdd_q)
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", year)] 
    
    #sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", (year))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}
  
  output3[[as.character(model)]] <- output
  
  
}

future_no_ac_cons <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons <- exp(future_no_ac_cons)
future_no_ac_cons[future_no_ac_cons<0] <- 0

future_impact_ac <- as.data.frame(future_ac_cons) - as.data.frame(future_no_ac_cons)
future_impact_ac[future_impact_ac<=0] <- NA

############################################

rowMedians <- function(data) {
  apply(data, 1, median, na.rm = TRUE)
}

# Group by ID and calculate average for col1 and col2 columns

future_acc <- data.frame(id=1:nrow(future_ac_adoption))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
col1_columns <- grepl(y, names(future_ac_adoption)) & grepl(x, names(future_ac_adoption))

varname <- paste0(x, ".", y)

future_acc <-bind_cols(future_acc, future_ac_adoption %>% dplyr::select(colnames(future_ac_adoption)[col1_columns]) %>% dplyr::summarise(!!varname := rowMedians(.)))

  }}

future_acc$id <- NULL

future_ac_adoption <- bind_cols(future_ac_adoption, future_acc)

####################################

# Group by ID and calculate average for col1 and col2 columns

future_acc <- data.frame(id=1:nrow(future_ac_cons))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(future_ac_cons)) & grepl(x, names(future_ac_cons))
    
    varname <- paste0(x, ".", y)
    
    future_acc <-bind_cols(future_acc, future_ac_cons %>% dplyr::select(colnames(future_ac_cons)[col1_columns]) %>% dplyr::summarise(!!varname := rowMedians(.)))
    
  }}

future_acc$id <- NULL

future_ac_cons <- bind_cols(future_ac_cons, future_acc)

####################################

# Group by ID and calculate average for col1 and col2 columns

future_acc <- data.frame(id=1:nrow(future_impact_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(future_impact_ac)) & grepl(x, names(future_impact_ac))
    
    varname <- paste0(x, ".", y)
    
    future_acc <-bind_cols(future_acc, future_impact_ac %>% dplyr::select(colnames(future_impact_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowMedians(.)))
    
  }}

future_acc$id <- NULL

future_impact_ac <- bind_cols(future_impact_ac, future_acc)


####################################

shape_ac <- bind_cols(shape, future_ac_adoption)

shape_ely_noac <- bind_cols(shape_all %>% st_set_geometry(NULL), future_no_ac_cons)

shape_ely_ac <- bind_cols(shape_all %>% st_set_geometry(NULL), future_impact_ac)

future_impact_ac_pctg <- future_impact_ac / as.data.frame(future_ac_cons)

shape_ely_diff <- bind_cols(shape_all %>% st_set_geometry(NULL), future_impact_ac)

shape_ely_pctg <- bind_cols(shape_all %>% st_set_geometry(NULL), future_impact_ac_pctg)

#######################

# produce total: multiply by grid cell population and AC penetration rate

shape_ely_diff$ely_total_SSP2_2010 <- (shape_ely_diff$SSP2.2010 * shape_ac$SSP2.2010 * (shape_ely_diff$pop_ssps_data_hist / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP2_2020 <- (shape_ely_diff$SSP2.2020 * shape_ac$SSP2.2020 * (shape_ely_diff$pop_ssps_data_ssp2_2020 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP2_2030 <- (shape_ely_diff$SSP2.2030 * shape_ac$SSP2.2030 * (shape_ely_diff$pop_ssps_data_ssp2_2030 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP2_2040 <- (shape_ely_diff$SSP2.2040 * shape_ac$SSP2.2040 * (shape_ely_diff$pop_ssps_data_ssp2_2040 / shape_ac$hhsize)) / 1e6

# shape_ely_diff$ely_total_SSP2_2020 <- ifelse(shape_ely_diff$ely_total_SSP2_2020 == 0, NA, shape_ely_diff$ely_total_SSP2_2020 )

shape_ely_diff$ely_total_SSP2_2050 <- (shape_ely_diff$SSP2.2050 * shape_ac$SSP2.2050 * (shape_ely_diff$pop_ssps_data_ssp2_2050 / shape_ac$hhsize)) / 1e6

# shape_ely_diff$ely_total_SSP5_2050 <- ifelse(shape_ely_diff$ely_total_SSP5_2050 == 0, NA, shape_ely_diff$ely_total_SSP5_2050 )

shape_ely_diff$ely_total_SSP5_2010 <- (shape_ely_diff$SSP5.2010 * shape_ac$SSP5.2010 * (shape_ely_diff$pop_ssps_data_hist / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP5_2020 <- (shape_ely_diff$SSP5.2020 * shape_ac$SSP5.2020 * (shape_ely_diff$pop_ssps_data_ssp5_2020 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP5_2030 <- (shape_ely_diff$SSP5.2030 * shape_ac$SSP5.2030 * (shape_ely_diff$pop_ssps_data_ssp5_2030 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP5_2040 <- (shape_ely_diff$SSP5.2040 * shape_ac$SSP5.2040 * (shape_ely_diff$pop_ssps_data_ssp5_2040 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP5_2050 <- (shape_ely_diff$SSP5.2050 * shape_ac$SSP5.2050 * (shape_ely_diff$pop_ssps_data_ssp5_2050 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP1_2010 <- (shape_ely_diff$SSP1.2010 * shape_ac$SSP1.2010 * (shape_ely_diff$pop_ssps_data_hist / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP1_2020 <- (shape_ely_diff$SSP1.2020 * shape_ac$SSP1.2020 * (shape_ely_diff$pop_ssps_data_ssp1_2020 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP1_2030 <- (shape_ely_diff$SSP1.2030 * shape_ac$SSP1.2030 * (shape_ely_diff$pop_ssps_data_ssp1_2030 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP1_2040 <- (shape_ely_diff$SSP1.2040 * shape_ac$SSP1.2040 * (shape_ely_diff$pop_ssps_data_ssp1_2040 / shape_ac$hhsize)) / 1e6

# shape_ely_diff$ely_total_SSP1_2020 <- ifelse(shape_ely_diff$ely_total_SSP1_2020 == 0, NA, shape_ely_diff$ely_total_SSP1_2020 )

shape_ely_diff$ely_total_SSP1_2050 <- (shape_ely_diff$SSP1.2050 * shape_ac$SSP1.2050 * (shape_ely_diff$pop_ssps_data_ssp1_2050 / shape_ac$hhsize)) / 1e6

# shape_ely_diff$ely_total_SSP3_2050 <- ifelse(shape_ely_diff$ely_total_SSP3_2050 == 0, NA, shape_ely_diff$ely_total_SSP3_2050 )

shape_ely_diff$ely_total_SSP3_2010 <- (shape_ely_diff$SSP3.2010 * shape_ac$SSP3.2010 * (shape_ely_diff$pop_ssps_data_hist / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP3_2020 <- (shape_ely_diff$SSP3.2020 * shape_ac$SSP3.2020 * (shape_ely_diff$pop_ssps_data_ssp3_2020 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP3_2030 <- (shape_ely_diff$SSP3.2030 * shape_ac$SSP3.2030 * (shape_ely_diff$pop_ssps_data_ssp3_2030 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP3_2040 <- (shape_ely_diff$SSP3.2040 * shape_ac$SSP3.2040 * (shape_ely_diff$pop_ssps_data_ssp3_2040 / shape_ac$hhsize)) / 1e6

shape_ely_diff$ely_total_SSP3_2050 <- (shape_ely_diff$SSP3.2050 * shape_ac$SSP3.2050 * (shape_ely_diff$pop_ssps_data_ssp3_2050 / shape_ac$hhsize)) / 1e6

# shape_ely_diff$ely_total_SSP5_2050 <- ifelse(shape_ely_diff$ely_total_SSP5_2050 == 0, NA, shape_ely_diff$ely_total_SSP5_2050 )

shape_ely_diff$geometry <- shape_all$geometry
shape_ely_diff <- st_as_sf(shape_ely_diff)

r0 <- fasterize(shape_ely_diff, r, "ely_total_SSP2_2020")
values(r0) <- ifelse(values(r0)<0, NA, values(r0))

r1 <- fasterize(shape_ely_diff, r, "ely_total_SSP2_2050")
values(r1) <- ifelse(values(r1)<0, NA, values(r1))

r2 <- fasterize(shape_ely_diff, r, "ely_total_SSP5_2050")
values(r2) <- ifelse(values(r2)<0, NA, values(r2))

r3 <- fasterize(shape_ely_diff, r, "ely_total_SSP1_2050")
values(r3) <- ifelse(values(r3)<0, NA, values(r3))

r4 <- fasterize(shape_ely_diff, r, "ely_total_SSP3_2050")
values(r4) <- ifelse(values(r4)<0, NA, values(r4))

#####################
#####################

# figures: SSPs2-5 at 2050 for both margins

newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

line = .8
cex = 1
side = 3
adj=-0.085


data("wrld_simpl")
wrld_simpl <- subset(wrld_simpl, wrld_simpl$ISO3!="ATA")


shape_ac$geometry <- shape_all$geometry
shape_ac <- st_as_sf(shape_ac)

shape_ely_diff$geometry <- shape_all$geometry
shape_ely_diff <- st_as_sf(shape_ely_diff)

###

pdf("results/graphs_tables/maps_ac.pdf", height = 3, width = 5*2.2)

par(mfrow = c(1, 3))

plot(fasterize(shape_ac, r, "SSP2.2020")*100, col = terrain.colors(length(seq(0, 1, by = .1)*100)-1, rev=T), main="AC penetration (% of households), 2020", breaks= seq(0, 1, by = .1)*100, ylim=range(-60:90) ,xaxt = "n", yaxt = "n")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("A", side=side, line=line, cex=cex, adj=adj)

plot(fasterize(shape_ac, r, "SSP2.2050")*100, col = terrain.colors(length(seq(0, 1, by = .1)*100)-1, rev=T), main="AC penetration (% of households), SSP245, 2050", breaks= seq(0, 1, by = .1)*100, ylim=range(-60:90),xaxt = "n", yaxt = "n")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("B", side=side, line=line, cex=cex, adj=adj)

plot(fasterize(shape_ac, r, "SSP5.2050")*100, col = terrain.colors(length(seq(0, 1, by = .1)*100)-1, rev=T), main="AC penetration (% of households), SSP585, 2050", breaks= seq(0, 1, by = .1)*100, ylim=range(-60:90) ,xaxt = "n", yaxt = "n")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("C", side=side, line=line, cex=cex, adj=adj)

dev.off()

###

pdf("results/graphs_tables/maps_ely.pdf", height = 3, width = 5*2.2)

par(mfrow = c(1, 3))

f1 = fasterize(shape_ely_diff, r, "ely_total_SSP2_2020")
f1_c = f1
f1_c[f1_c<0.1] <- NA
f1_c[f1 < 1] <- 1
f1_c[f1 > 1 & f1 < 10] <- 2
f1_c[f1 > 10 & f1 < 100] <- 3
f1_c[f1 > 100 & f1 < 1000] <- 4
f1_c[f1 > 1000] <- 5
f1_c <- ratify(f1_c) # https://stackoverflow.com/questions/23840178/how-to-write-a-raster-with-rat-factors-in-r-raster-package
rat <- levels(f1_c)[[1]]#get the values of the unique cell frot the attribute table
rat$legend <- c('<1', "1-10", "10-100", "100-1000", "1000+")[1:nrow(levels(f1_c)[[1]])]
levels(f1_c) <- rat

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=heat.colors(5, rev=T), main="AC electricity use (GWh/yr.), 2010", ylim=range(-60:90) ,xaxt = "n", yaxt = "n", legend=F)
#legend(x='bottom', legend =rat$legend,fill = heat.colors(5, rev=T), horiz=T, box.col="transparent")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("D", side=side, line=line, cex=cex, adj=adj)

f1 = fasterize(shape_ely_diff, r, "ely_total_SSP2_2050")
f1_c = f1
f1_c[f1_c<0.1] <- NA
f1_c[f1 < 1] <- 1
f1_c[f1 > 1 & f1 < 10] <- 2
f1_c[f1 > 10 & f1 < 100] <- 3
f1_c[f1 > 100 & f1 < 1000] <- 4
f1_c[f1 > 1000] <- 5
f1_c <- ratify(f1_c) # https://stackoverflow.com/questions/23840178/how-to-write-a-raster-with-rat-factors-in-r-raster-package
rat <- levels(f1_c)[[1]]#get the values of the unique cell frot the attribute table
rat$legend <- c('<1', "1-10", "10-100", "100-1000", "1000+")[1:nrow(levels(f1_c)[[1]])]
levels(f1_c) <- rat

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=heat.colors(5, rev=T), main="AC electricity use (GWh/yr.), 2050, SSP245", ylim=range(-60:90) ,xaxt = "n", yaxt = "n", legend=F)
legend(x='bottom', legend =rat$legend,fill = heat.colors(5, rev=T), horiz=T, box.col="transparent")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("E", side=side, line=line, cex=cex, adj=adj)


f1 = fasterize(shape_ely_diff, r, "ely_total_SSP5_2050")
f1_c = f1
f1_c[f1_c<0.1] <- NA
f1_c[f1 < 1] <- 1
f1_c[f1 > 1 & f1 < 10] <- 2
f1_c[f1 > 10 & f1 < 100] <- 3
f1_c[f1 > 100 & f1 < 1000] <- 4
f1_c[f1 > 1000] <- 5
f1_c <- ratify(f1_c) # https://stackoverflow.com/questions/23840178/how-to-write-a-raster-with-rat-factors-in-r-raster-package
rat <- levels(f1_c)[[1]]#get the values of the unique cell frot the attribute table
rat$legend <- c('<1', "1-10", "10-100", "100-1000", "1000+")[1:nrow(levels(f1_c)[[1]])]
levels(f1_c) <- rat

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=heat.colors(5, rev=T), main="AC electricity use (GWh/yr.), 2050, SSP585", ylim=range(-60:90) ,xaxt = "n", yaxt = "n", legend=F)
#legend(x='bottom', legend =rat$legend,fill = heat.colors(5, rev=T), horiz=T, box.col="transparent")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("F", side=side, line=line, cex=cex, adj=adj)

dev.off()


###################
###################

pop_ssps <- list.files(path=paste0(wd, "/supporting_data/pop_downscaled_spps"), recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- stack(pop_ssps[2])[[10]]

#

source("maps_gridded_v3.R")


#########################

write_rds(shape_ac, file = "output_data/shape_ac.Rds")
write_rds(shape_ely_diff, file = "output_data/shape_ely_diff.Rds")

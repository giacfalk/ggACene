### 1 (baseline)

model <- cmip6_models[1]

output_baseline <- list()

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  orig_data_bk <- shape
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    #orig_data$weight = 0
    
    orig_data$ln_total_exp_usd_2011 = shape[,paste0("GDP_", ssp, "_", "2020")]
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output_baseline[[as.character(ssp)]] <- output2
  
}

future_ac_adoption_baseline <- do.call(cbind, lapply(output_baseline, as.data.frame))
future_ac_adoption_baseline[future_ac_adoption_baseline<=0] <- NA
shape_ac_baseline <- bind_cols(shape, future_ac_adoption_baseline)

##################
##################

### 2 (climate only)

output_cc <- list()

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  orig_data_bk <- shape
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    #orig_data$weight = 0
    
    orig_data$ln_total_exp_usd_2011 = shape[,paste0("GDP_", ssp, "_", "2020")]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
    
  }
  
  output_cc[[as.character(ssp)]] <- output2
  
}

future_ac_adoption_cc <- do.call(cbind, lapply(output_cc, as.data.frame))
future_ac_adoption_cc[future_ac_adoption_cc<=0] <- NA
shape_ac_cc <- bind_cols(shape, future_ac_adoption_cc)

##################
##################

### 3 (econ growth only)

output_econgrowth <- list()

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
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
     
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
    
  }
  
  output_econgrowth[[as.character(ssp)]] <- output2
  
}

future_ac_adoption_econgrowth <- do.call(cbind, lapply(output_econgrowth, as.data.frame))
future_ac_adoption_econgrowth[future_ac_adoption_econgrowth<=0] <- NA
shape_ac_econgrowth <- bind_cols(shape, future_ac_adoption_econgrowth)

##################
##################

### 4 (social drivers only)

output_socialdrivers <- list()

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  orig_data_bk <- shape
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    #orig_data$weight = 0
    
    orig_data$ln_total_exp_usd_2011 = shape[,paste0("GDP_", ssp, "_", "2020")]
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)]
    
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
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", year)]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
    
  }
  
  output_socialdrivers[[as.character(ssp)]] <- output2
  
}

future_ac_adoption_socialdrivers <- do.call(cbind, lapply(output_socialdrivers, as.data.frame))
future_ac_adoption_socialdrivers[future_ac_adoption_socialdrivers<=0] <- NA
shape_ac_socialdrivers <- bind_cols(shape, future_ac_adoption_socialdrivers)

##################
##################

### 5 (population only)

output_population <- list()

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  orig_data_bk <- shape
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    #orig_data$weight = 0
    
    orig_data$ln_total_exp_usd_2011 = shape[,paste0("GDP_", ssp, "_", "2020")]
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, "2020", "2020"), "_", model)]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", "2020")]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ac_model, orig_data, type="prob", na.action = na.omit)[,2]
    
    output2[[as.character(year)]] <- projected
    
    
  }
  
  output_population[[as.character(ssp)]] <- output2
  
}

future_ac_adoption_population <- do.call(cbind, lapply(output_population, as.data.frame))
future_ac_adoption_population[future_ac_adoption_population<=0] <- NA
shape_ac_population <- bind_cols(shape, future_ac_adoption_population)

###################################
###################################

save(shape, shape_ac_baseline, shape_ac_cc, shape_ac_econgrowth, shape_ac_socialdrivers, shape_ac_population, file="fordriver_decompo.Rdata")

load("fordriver_decompo.Rdata")

# plot

shape_ac_baseline <- dplyr::select(shape_ac_baseline, id, 4201:4220)
shape_ac_cc <- dplyr::select(shape_ac_cc, id,  4201:4220)
shape_ac_econgrowth <- dplyr::select(shape_ac_econgrowth, id,  4201:4220)
shape_ac_socialdrivers <- dplyr::select(shape_ac_socialdrivers, id,  4201:4220)
shape_ac_population <- dplyr::select(shape_ac_population, id,  4201:4220)

pops <- dplyr::select(shape, id,  starts_with("pop_ssps_data"))  %>% dplyr::select(!contains("ssp4"))
pops$pop_ssps_data_hist <- NULL
pops <- melt(pops, "id")
colnames(pops)[3] <- "population"

pops$variable <- gsub("pop_ssps_data_", "", pops$variable)
pops$variable <- gsub("_", ".X", pops$variable)
pops$variable <- toupper(pops$variable)

shape_ac_baseline <- melt(shape_ac_baseline, "id")
shape_ac_baseline <- merge(shape_ac_baseline, pops, by=c("id", "variable"))

shape_ac_cc <- melt(shape_ac_cc, "id")
shape_ac_cc <- merge(shape_ac_cc, pops, by=c("id", "variable"))

shape_ac_econgrowth <- melt(shape_ac_econgrowth, "id")
shape_ac_econgrowth <- merge(shape_ac_econgrowth, pops, by=c("id", "variable"))

shape_ac_socialdrivers <- melt(shape_ac_socialdrivers, "id")
shape_ac_socialdrivers <- merge(shape_ac_socialdrivers, pops, by=c("id", "variable"))

shape_ac_population <- melt(shape_ac_population, "id")
shape_ac_population <- merge(shape_ac_population, pops, by=c("id", "variable"))

#

shape_ac_baseline <- shape_ac_baseline %>% group_by(variable) %>% dplyr::summarise(value=weighted.mean(value, population, na.rm=T))
shape_ac_cc <- shape_ac_cc %>% group_by(variable) %>% dplyr::summarise(value=weighted.mean(value, population, na.rm=T))
shape_ac_econgrowth <- shape_ac_econgrowth %>% group_by(variable) %>% dplyr::summarise(value=weighted.mean(value, population, na.rm=T))
shape_ac_socialdrivers <- shape_ac_socialdrivers %>% group_by(variable) %>% dplyr::summarise(value=weighted.mean(value, population, na.rm=T))
shape_ac_population <- shape_ac_population %>% group_by(variable) %>% dplyr::summarise(value=weighted.mean(value, population, na.rm=T))

shape_ac_baseline$driver <- " historical"
shape_ac_cc$driver <- "cc"
shape_ac_econgrowth$driver <- "econgrowth"
shape_ac_socialdrivers$driver <- "socialdrivers"
shape_ac_population$driver <- "population"

shape_ac_cc$value <- shape_ac_cc$value - shape_ac_baseline$value
shape_ac_econgrowth$value <- shape_ac_econgrowth$value - shape_ac_baseline$value
shape_ac_socialdrivers$value <- shape_ac_socialdrivers$value - shape_ac_baseline$value
shape_ac_population$value <- shape_ac_population$value - shape_ac_baseline$value

###

ss <- bind_rows(shape_ac_cc, shape_ac_econgrowth, shape_ac_socialdrivers, shape_ac_baseline)

ss$scenario <- substr(ss$variable, 1, 4)
ss$year <- substr(ss$variable, 7, 10)

ss$driver <- factor(ss$driver, levels=rev(c(" historical", "cc", "econgrowth", "socialdrivers")))

deco_ac <- ggplot(ss %>% dplyr::filter(year>2010))+
  geom_col(aes(x=year, y=value*100, fill=driver))+
  facet_wrap(vars(scenario))+
  ylab("Global AC penetration rate (%)")+
  theme_classic()+
  xlab("")+
  scale_fill_brewer(palette = "Set3")


ggsave("results/graphs_tables/deco_ac.png", deco_ac, scale=1.3, height = 4, width =4)


###################################
###################################
###################################
###################################

model <- cmip6_models[1]

load("fordriver_decompo.Rdata")

### 1 (baseline)

shape_ac_baseline <- na_mean(shape_ac_baseline)

orig_data <- shape_ac_baseline
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

# loop for all ssps and time-steps
  
  for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
    
    print(ssp)
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
    
    output2 <- list()
    
    for (year in seq(2010, 2050, 10)){
      
      orig_data <- orig_data_bk
      
      print(year)
      
      orig_data$phat0_obs = 1
      
      orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
      
      orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
      
      orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
      
      orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
      
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
      orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
      
      #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
      
      orig_data$macroregion  <- shape$macroregion 
      
      orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
      
      orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
      
      if (year == 2010){
        
        orig_data$weight  = shape$pop_ssps_data_hist
        
        
      } else{
        
        orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
        
      }
      
      #
      projected <- predict(ely_model, orig_data)
      
      output2[[as.character(year)]] <- projected
      
    }
    
    output3[[as.character(ssp)]] <- output2
    
}

future_ac_cons_baseline <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons_baseline <- exp(future_ac_cons_baseline)
future_ac_cons_baseline[future_ac_cons_baseline<=0] <- NA

###########################

output3 <- list()

# loop for all ssps and time-steps

  for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
    
    print(ssp)
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
    
    output2 <- list()
    
    for (year in seq(2010, 2050, 10)){
      
      orig_data <- orig_data_bk
      
      print(year)
      
      orig_data$phat0_obs = 0
      
      orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
      
      orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
      
      orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
      
      orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
      
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
      
      orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
      
      #sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
      
      orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
      
      orig_data$macroregion  <- shape$macroregion 
      
      orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
      
      if (year == 2010){
        
        orig_data$weight  = shape$pop_ssps_data_hist
        
        
      } else{
        
        orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
        
      }
      
      #
      projected <- predict(ely_model, orig_data)
      
      output2[[as.character(year)]] <- projected
      
    }
    
    output3[[as.character(ssp)]] <- output2
  
}

future_no_ac_cons_baseline <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons_baseline <- exp(future_no_ac_cons_baseline)
future_no_ac_cons_baseline[future_no_ac_cons_baseline<0] <- 0

future_impact_ac_baseline <- as.data.frame(future_ac_cons_baseline) - as.data.frame(future_no_ac_cons_baseline)
future_impact_ac_baseline[future_impact_ac_baseline<=0] <- NA
shape_ely_baseline <- bind_cols(shape, future_impact_ac_baseline)

##############

### 2 (climate only)

shape_ac_cc <- na_mean(shape_ac_cc)

orig_data <- shape_ac_cc
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 1 
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
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
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_ac_cons_cc <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons_cc <- exp(future_ac_cons_cc)
future_ac_cons_cc[future_ac_cons_cc<=0] <- NA

###########################

output3 <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 0
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_no_ac_cons_cc <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons_cc <- exp(future_no_ac_cons_cc)
future_no_ac_cons_cc[future_no_ac_cons_cc<0] <- 0

future_impact_ac_cc <- as.data.frame(future_ac_cons_cc) - as.data.frame(future_no_ac_cons_cc)
future_impact_ac_cc[future_impact_ac_cc<=0] <- NA
shape_ely_cc <- bind_cols(shape, future_impact_ac_cc)

##############

### 2 (economic growth only)

shape_ac_econgrowth <- na_mean(shape_ac_econgrowth)

orig_data <- shape_ac_econgrowth
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 1
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", year)])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_ac_cons_econgrowth <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons_econgrowth <- exp(future_ac_cons_econgrowth)
future_ac_cons_econgrowth[future_ac_cons_econgrowth<=0] <- NA

###########################

output3 <- list()

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
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_no_ac_cons_econgrowth <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons_econgrowth <- exp(future_no_ac_cons_econgrowth)
future_no_ac_cons_econgrowth[future_no_ac_cons_econgrowth<0] <- 0

future_impact_ac_econgrowth <- as.data.frame(future_ac_cons_econgrowth) - as.data.frame(future_no_ac_cons_econgrowth)
future_impact_ac_econgrowth[future_impact_ac_econgrowth<=0] <- NA
shape_ely_econgrowth <- bind_cols(shape, future_impact_ac_econgrowth)

##############

### 3 (social drivers only)

shape_ac_socialdrivers <- na_mean(shape_ac_socialdrivers)

orig_data <- shape_ac_socialdrivers
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 1
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_ac_cons_socialdrivers <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons_socialdrivers <- exp(future_ac_cons_socialdrivers)
future_ac_cons_socialdrivers[future_ac_cons_socialdrivers<=0] <- NA

###########################

output3 <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 0
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", ("2020"))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_no_ac_cons_socialdrivers <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons_socialdrivers <- exp(future_no_ac_cons_socialdrivers)
future_no_ac_cons_socialdrivers[future_no_ac_cons_socialdrivers<0] <- 0

future_impact_ac_socialdrivers <- as.data.frame(future_ac_cons_socialdrivers) - as.data.frame(future_no_ac_cons_socialdrivers)
future_impact_ac_socialdrivers[future_impact_ac_socialdrivers<=0] <- NA
shape_ely_socialdrivers <- bind_cols(shape, future_impact_ac_socialdrivers)

### 4 (population only)

shape_ac_population <- na_mean(shape_ac_population)

orig_data <- shape_ac_population
orig_data$geometry <- NULL
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

global$ac = ifelse(global$ac=="Yes", 1, 0)

###

output3 <- list()

# loop for all ssps and time-steps

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 1
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #orig_data$sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_ac_cons_population <- do.call(cbind, lapply(output3, as.data.frame))
future_ac_cons_population <- exp(future_ac_cons_population)
future_ac_cons_population[future_ac_cons_population<=0] <- NA

###########################

output3 <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", ifelse(ssp=="SSP1", "rcp26", ifelse(ssp=="SSP3", "rcp70", "rcp85")))
  
  output2 <- list()
  
  for (year in seq(2010, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$phat0_obs = 0
    
    orig_data$ln_total_exp_usd_2011 = (shape[,paste0("GDP_", ssp, "_", "2020")])
    
    orig_data$mean_CDD18_db  = log(shape[,paste0("CDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$mean_HDD18_db  = log(shape[,paste0("HDD_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2015, "2020"), "_", model)] +1)
    
    orig_data$hurs  = shape[,paste0("hurs_", ifelse(ssp=="SSP2", 245, ifelse(ssp=="SSP1", 126, ifelse(ssp=="SSP3", 370, 585))), "_", ifelse(year==2010, 2020, "2020"), "_", model)]
    
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
    
    orig_data$edu_head_2 <- shape[,paste0("edu_", ssp, "_", "2020")] 
    
    #sex_head  <- shape[,paste0("women_", ssp, "_", year)] 
    
    orig_data$age_head <- shape[,paste0("age_", ssp, "_", "2020")] 
    
    orig_data$macroregion  <- shape$macroregion 
    
    orig_data$urban_sh  = shape[,paste0("URB_", ssp, "_", ("2020"))]
    
    if (year == 2010){
      
      orig_data$weight  = shape$pop_ssps_data_hist
      
      
    } else{
      
      orig_data$weight  = shape[,paste0("pop_ssps_data_", tolower(ssp), "_", (year))]
      
    }
    
    #
    projected <- predict(ely_model, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output3[[as.character(ssp)]] <- output2
  
}

future_no_ac_cons_population <- do.call(cbind, lapply(output3, as.data.frame))
future_no_ac_cons_population <- exp(future_no_ac_cons_population)
future_no_ac_cons_population[future_no_ac_cons_population<0] <- 0

future_impact_ac_population <- as.data.frame(future_ac_cons_population) - as.data.frame(future_no_ac_cons_population)
future_impact_ac_population[future_impact_ac_population<=0] <- NA
shape_ely_population <- bind_cols(shape, future_impact_ac_population)

###############################
###############################

save(shape, shape_ely_baseline, shape_ely_cc, shape_ely_econgrowth, shape_ely_socialdrivers, shape_ely_population, file="fordriver_decompo_ely.Rdata")

# plot

load("fordriver_decompo_ely.Rdata")

# plot

shape_ely_baseline <- dplyr::select(shape_ely_baseline, id, 4201:4220)
shape_ely_cc <- dplyr::select(shape_ely_cc, id,  4201:4220)
shape_ely_econgrowth <- dplyr::select(shape_ely_econgrowth, id,  4201:4220)
shape_ely_socialdrivers <- dplyr::select(shape_ely_socialdrivers, id,  4201:4220)
shape_ely_population <- dplyr::select(shape_ely_population, id,  4201:4220)

pops <- dplyr::select(shape, id,  starts_with("pop_ssps_data"))  %>% dplyr::select(!contains("ssp4"))
pops$pop_ssps_data_hist <- NULL
pops <- melt(pops, "id")
colnames(pops)[3] <- "population"

hhsize <- dplyr::select(shape, id,  hhsize)
hhsize <- melt(hhsize, "id")
colnames(hhsize)[3] <- "hhsize"
hhsize$variable <- NULL

pops$variable <- gsub("pop_ssps_data_", "", pops$variable)
pops$variable <- gsub("_", ".X", pops$variable)
pops$variable <- toupper(pops$variable)

shape_ely_baseline <- melt(shape_ely_baseline, "id")
shape_ely_baseline <- merge(shape_ely_baseline, pops, by=c("id", "variable"))
shape_ely_baseline <- merge(shape_ely_baseline, hhsize, by=c("id"))

shape_ely_cc <- melt(shape_ely_cc, "id")
shape_ely_cc <- merge(shape_ely_cc, pops, by=c("id", "variable"))
shape_ely_cc <- merge(shape_ely_cc, hhsize, by=c("id"))

shape_ely_econgrowth <- melt(shape_ely_econgrowth, "id")
shape_ely_econgrowth <- merge(shape_ely_econgrowth, pops, by=c("id", "variable"))
shape_ely_econgrowth <- merge(shape_ely_econgrowth, hhsize, by=c("id"))

shape_ely_socialdrivers <- melt(shape_ely_socialdrivers, "id")
shape_ely_socialdrivers <- merge(shape_ely_socialdrivers, pops, by=c("id", "variable"))
shape_ely_socialdrivers <- merge(shape_ely_socialdrivers, hhsize, by=c("id"))

shape_ely_population <- melt(shape_ely_population, "id")
shape_ely_population <- merge(shape_ely_population, pops, by=c("id", "variable"))
shape_ely_population <- merge(shape_ely_population, hhsize, by=c("id"))

#

weighted_sum <- function(values, weights) {
  return(sum(values * weights, na.rm=T))
}

#

shape_ely_baseline <- shape_ely_baseline %>% group_by(variable) %>% dplyr::summarise(value=weighted_sum(value, population/hhsize))
shape_ely_cc <- shape_ely_cc%>% group_by(variable) %>% dplyr::summarise(value=weighted_sum(value, population/hhsize))
shape_ely_econgrowth <- shape_ely_econgrowth %>% group_by(variable) %>% dplyr::summarise(value=weighted_sum(value, population/hhsize))
shape_ely_socialdrivers <- shape_ely_socialdrivers %>% group_by(variable) %>% dplyr::summarise(value=weighted_sum(value, population/hhsize))
shape_ely_population <- shape_ely_population %>% group_by(variable) %>% dplyr::summarise(value=weighted_sum(value, population/hhsize))

shape_ely_baseline$driver <- " historical"
shape_ely_cc$driver <- "cc"
shape_ely_econgrowth$driver <- "econgrowth"
shape_ely_socialdrivers$driver <- "socialdrivers"
shape_ely_population$driver <- "population"

shape_ely_cc$value <- shape_ely_cc$value - shape_ely_baseline$value
shape_ely_econgrowth$value <- shape_ely_econgrowth$value - shape_ely_baseline$value
shape_ely_socialdrivers$value <- shape_ely_socialdrivers$value - shape_ely_baseline$value
shape_ely_population$value <- shape_ely_population$value - shape_ely_baseline$value

###

ss <- bind_rows(shape_ely_cc, shape_ely_econgrowth, shape_ely_socialdrivers, shape_ely_baseline)

ss$scenario <- substr(ss$variable, 1, 4)
ss$year <- substr(ss$variable, 7, 10)

ss$driver <- factor(ss$driver, levels=rev(c(" historical", "cc", "econgrowth", "socialdrivers")))

ss$value <- ifelse(ss$value<0, -ss$value, ss$value)

ss$value[ss$driver==" historical"] <- 0.55e+12

deco_ely <- ggplot(ss %>% dplyr::filter(year>2010))+
  geom_col(aes(x=year, y=value/1e9, fill=driver))+
  facet_wrap(vars(scenario))+
  ylab("Global AC electricity use (TWh/yr.)")+
  theme_classic()+
  xlab("")+
  scale_fill_brewer(palette = "Set3")


ggsave("results/graphs_tables/deco_ely.png", deco_ely, scale=1.3, height = 4, width =4)

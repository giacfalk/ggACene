## 1) Load libraries and data ##
library(tidyverse)
library(reshape2)
library(sf)
library(gg.layers)
####

setwd(wd)

load(paste0(wd, "supporting_data/data_for_global_spline_v2.Rds"))

############

shape <- st_as_sf(shape)
shape <- filter(shape, !is.na(region))

shape_gdp <- dplyr::select(shape, contains("GDP"), region)
shape_gdp$geometry <- NULL

shape_gdp <- melt(shape_gdp, 26)
shape_gdp$ssp <- substr(shape_gdp$variable, 5, 8)
shape_gdp$year <- substr(shape_gdp$variable, 10, 13)

shape_gdp <- filter(shape_gdp, year=="2020" | year=="2050")
shape_gdp <- filter(shape_gdp, exp(value)<100000)

shape_gdp <- filter(shape_gdp, ssp!="SSP4")

ggplot(shape_gdp)+
  theme_classic()+
  gg.layers::geom_boxplot2(aes(y=exp(value),  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
    facet_wrap(vars(region), ncol=1)+
  ggsci::scale_fill_npg(name="Scenario and year")+
  ylab("PPP per-capita expenditure, 2011 USDs")+
  xlab("Scenario")+
  theme(legend.position = "none")

ggsave("expenditure_evolution.png", scale=2, width = 3, height = 5)

#

shape_2 <- shape
shape_2$region <- as.character(shape_2$region)
shape_2$region <- " GLOBAL"

shape$region <- as.character(shape$region)
shape <- bind_rows(shape, shape_2)
shape$region <- as.factor(shape$region)

shape_pop <- dplyr::select(shape, contains("pop"), region)
shape_pop$POP2005 = NULL
shape_pop$geometry <- NULL

shape_pop %>% group_by(region) %>% dplyr::summarise(pop_ssps_data_ssp2_2020=sum(pop_ssps_data_ssp2_2020)/1e9)

shape_pop <- melt(shape_pop, 27)
shape_pop$ssp <- gsub("pop_ssps_data_", "", shape_pop$variable)
shape_pop$ssp <- gsub("_2010|_2020|_2030|_2040|_2050", "", shape_pop$ssp)
shape_pop <- filter(shape_pop, ssp!="hist")

shape_pop$year <- gsub("pop_ssps_data_ssp1_", "", shape_pop$variable)
shape_pop$year <- gsub("pop_ssps_data_ssp4_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp3_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp2_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_ssp5_", "", shape_pop$year)
shape_pop$year <- gsub("pop_ssps_data_", "", shape_pop$year)

shape_pop <- filter(shape_pop, year=="2020" | year=="2050")

shape_pop <- filter(shape_pop, ssp!="ssp4")

ggplot(shape_pop %>% group_by(ssp, year, region) %>% dplyr::summarise(value=sum(value, na.rm=T)))+
  theme_classic()+
  geom_col(aes(y=value/1e9,  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
  facet_wrap(vars(region), ncol=1, scales = "free_y")+
  ggsci::scale_fill_npg(name="Scenario and year", na.translate=F)+
  ylab("Population count, billion")+
  xlab("Scenario")+
  theme(legend.position = "none")

ggsave("population_evolution.png", scale=2, width = 3, height = 5)

#


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

ggplot(shape_cdd)+
  theme_classic()+
  geom_boxplot2(aes(y=value,  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
  facet_wrap(vars(region), ncol=4, scales = "free")+
  ggsci::scale_fill_npg(name="Scenario and year")+
  ylab("Cooling Degree Days / yr")+
  theme(legend.position = "none")


ggsave("cdd_evolution.png", scale=2.5, width = 8, height = 3)

#

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

ggplot(shape_hdd)+
  theme_classic()+
  geom_boxplot2(aes(y=value,  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
  facet_wrap(vars(region), ncol=4, scales = "free")+
  ggsci::scale_fill_npg(name="Scenario and year")+
  ylab("Heating Degree Days / yr")+
  theme(legend.position = "none")


ggsave("hdd_evolution.png", scale=2.5, width = 8, height = 3)

#

shape_age <- dplyr::select(shape, contains("age"), region)
shape_age$geometry <- NULL

shape_age <- melt(shape_age, 26)
shape_age$ssp <- substr(shape_age$variable, 5, 8)
shape_age$year <- substr(shape_age$variable, 10, 13)

shape_age <- filter(shape_age, ssp!="SSP4")

shape_age <- filter(shape_age, year=="2020" | year=="2050")

shape_age <- filter(shape_age, !is.na(region))

shape_age <- shape_age %>% 
  group_by(region, ssp, year) %>%
  mutate(med_age = median(as.numeric(value), na.rm = T))

ggplot(shape_age)+
  theme_classic()+
  geom_boxplot2(aes(y=value,  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
  facet_wrap(vars(region), ncol=2, scales = "free")+
  ggsci::scale_fill_npg(name="Scenario and year")+
  ylab("Count of grid cells")+
  xlab("Age of population")+
  theme(legend.position = "none")

ggsave("age_evolution.png", scale=2, width = 8, height = 6)

#

shape_edu <- dplyr::select(shape, contains("edu"), region)
shape_edu$geometry <- NULL

shape_edu <- melt(shape_edu, 26)
shape_edu$ssp <- substr(shape_edu$variable, 5, 8)
shape_edu$year <- substr(shape_edu$variable, 10, 13)

shape_edu <- filter(shape_edu, ssp!="SSP4")

shape_edu <- filter(shape_edu, year=="2020" | year=="2050")

shape_edu <- filter(shape_edu, !is.na(region))

shape_edu <- shape_edu %>% 
  group_by(region, ssp, year) %>%
  mutate(med_edu = median(as.numeric(value), na.rm = T))

ggplot(shape_edu)+
  theme_classic()+
  geom_boxplot2(aes(y=value,  x=as.factor(interaction(ssp, year)), fill=as.factor(interaction(ssp, year))))+
  facet_wrap(vars(region), ncol=2, scales = "free")+
  ggsci::scale_fill_npg(name="Scenario and year")+
  ylab("Count of grid cells")+
  xlab("Education level of population")+
  theme(legend.position = "none")

ggsave("edu_evolution.png", scale=2, width = 8, height = 6)

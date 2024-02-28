rm(list=ls(all=TRUE)) # Removes all previously created variables
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

####

# Set users
user <- 'fp'
user <- 'gf'
user <- 'gf_server'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

if (user=='gf_server') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

setwd(paste0(stub, "rscripts/global_spline"))

##################################################

data <- paste(stub,'results/regressions/for_projections', sep='')
output <- paste(stub,'results/regressions/', sep='')

# Load country data
load("results/global_wgt_dmcf.RData")

global = reg_ely$data

global <- as.data.frame(global)

###

l <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/climate/processed/hurs", pattern="rds", full.names = T)
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

global <- dplyr::select(global, country, ac, ln_ely_q, country, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban_sh, edu_head_2,  age_head, n_members, weight, hurs, ely_p_usd_2011)

global <- dplyr::select(global, where(not_all_na))

gc()

######################

# adjustment expenditure to gdp per capita

adj_gtap <- read.csv("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/global_spline/adj_gtap/exp_gdp_capita_markups.csv")

adj_gtap$REG <- toupper(adj_gtap$REG)

global$iso3c <- countrycode::countrycode(global$country, 'country.name', 'iso3c')

global <- merge(global, adj_gtap, by.x="iso3c", by.y="REG")

#global$ln_total_exp_usd_2011 <- log(exp(global$ln_total_exp_usd_2011) / global$markup)

# add a macroregion variable

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

load(paste0(stub, "6-Projections/rscripts/global_spline/supporting_data/adj_factors.Rds"))

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

##########
#########

# calibrate gridded input data

keep(stub, global, sure=T)

load("supporting_data/data_for_global_spline_v2.Rds")


# summary(global %>% filter(country=="France") %>% mutate(ln_total_exp_usd_2011 = (exp(ln_total_exp_usd_2011)*(weight))/mean(weight, na.rm=T)) %>% dplyr::select(ln_total_exp_usd_2011) %>% pull())
# 
# summary(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010))
# 
# summary(ifelse(((as.numeric(scale(na.omit(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010)), center = median(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010), na.rm=T)))) + 1) *  median(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010), na.rm=T)<0, 1, (as.numeric(scale(na.omit(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010)), center = median(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010), na.rm=T)))) + 1) *  median(exp(filter(shape, ISO3=="FRA")$GDP_SSP2_2010), na.rm=T))
# 
# 
# ####
# 
# sm <- ifelse(((as.numeric(scale(na.omit(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010)), center = median(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010), na.rm=T)))) + 1) *  median(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010), na.rm=T)<0, 1, (as.numeric(scale(na.omit(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010)), center = median(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010), na.rm=T)))) + 1) *  median(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010), na.rm=T)
# 
# 
# ggplot()+
#   geom_density(aes(x=global %>% filter(country=="Germany") %>% mutate(ln_total_exp_usd_2011 = (exp(ln_total_exp_usd_2011)*(weight))/mean(weight, na.rm=T)) %>% dplyr::select(ln_total_exp_usd_2011) %>% pull()), colour="red")+
#   geom_density(aes(x=(exp(filter(shape, ISO3=="DEU")$GDP_SSP2_2010))), colour="blue")+
#   geom_density(aes(x=sm))+
#   coord_cartesian(xlim=c(0, 100000))


###############################

shape <- shape %>% group_by(ISO3) %>% mutate_at(vars(starts_with("GDP_")), list(~log(ifelse(((as.numeric(scale(exp(.), center = median(exp(.), na.rm=T)))) + 1) *  median(exp(.), na.rm=T)<0, 1, (as.numeric(scale(exp(.), center = median(exp(.), na.rm=T)))) + 1) *  median(exp(.), na.rm=T))))

#summary(exp(shape$GDP_SSP2_2010[shape$ISO3=="USA"]))

###
# calibrate total pop

# 

c <- read.csv("supporting_data/SspDb_compare_regions_2013-06-12.csv")
c <- filter(c, VARIABLE=="Population" & REGION=="World" & MODEL=="IIASA-WiC POP")
c$SCENARIO <- substr(c$SCENARIO, 1, 4)

totpop_pop_current <- c$X2010[c$SCENARIO=="SSP2"] * 1e6
totpop_pop_ssp2_2020 <-  c$X2020[c$SCENARIO=="SSP2"] * 1e6
totpop_pop_ssp2_2030 <-  c$X2030[c$SCENARIO=="SSP2"] * 1e6
totpop_pop_ssp2_2040 <-  c$X2040[c$SCENARIO=="SSP2"] * 1e6
totpop_pop_ssp2_2050 <-  c$X2050[c$SCENARIO=="SSP2"] * 1e6
totpop_pop_ssp5_2020 <-  c$X2020[c$SCENARIO=="SSP5"] * 1e6
totpop_pop_ssp5_2030 <-  c$X2030[c$SCENARIO=="SSP5"] * 1e6
totpop_pop_ssp5_2040 <-  c$X2040[c$SCENARIO=="SSP5"] * 1e6
totpop_pop_ssp5_2050 <-  c$X2050[c$SCENARIO=="SSP5"] * 1e6


#############
#############
#############

adjfact_pop_current <- totpop_pop_current / sum(shape$pop_ssps_data_hist, na.rm=T)
adjfact_pop_ssp2_2020 <-  totpop_pop_ssp2_2020 / sum(shape$pop_ssps_data_ssp2_2020, na.rm=T)
adjfact_pop_ssp2_2030 <-  totpop_pop_ssp2_2030 / sum(shape$pop_ssps_data_ssp2_2030, na.rm=T)
adjfact_pop_ssp2_2040 <-  totpop_pop_ssp2_2040 / sum(shape$pop_ssps_data_ssp2_2040, na.rm=T)
adjfact_pop_ssp2_2050 <-  totpop_pop_ssp2_2050 / sum(shape$pop_ssps_data_ssp2_2050, na.rm=T)
adjfact_pop_ssp5_2020 <-  totpop_pop_ssp5_2020 / sum(shape$pop_ssps_data_ssp5_2020, na.rm=T)
adjfact_pop_ssp5_2030 <-  totpop_pop_ssp5_2030 / sum(shape$pop_ssps_data_ssp5_2030, na.rm=T)
adjfact_pop_ssp5_2040 <-  totpop_pop_ssp5_2040 / sum(shape$pop_ssps_data_ssp5_2040, na.rm=T)
adjfact_pop_ssp5_2050 <-  totpop_pop_ssp5_2050 / sum(shape$pop_ssps_data_ssp5_2050, na.rm=T)

#############
#############
#############

shape$pop_ssps_data_hist <- shape$pop_ssps_data_hist * adjfact_pop_current
shape$pop_ssps_data_ssp2_2020 <- shape$pop_ssps_data_ssp2_2020 * adjfact_pop_ssp2_2020
shape$pop_ssps_data_ssp2_2030 <- shape$pop_ssps_data_ssp2_2030 * adjfact_pop_ssp2_2030
shape$pop_ssps_data_ssp2_2040 <- shape$pop_ssps_data_ssp2_2040 * adjfact_pop_ssp2_2040
shape$pop_ssps_data_ssp2_2050 <- shape$pop_ssps_data_ssp2_2050 * adjfact_pop_ssp2_2050
shape$pop_ssps_data_ssp5_2020 <- shape$pop_ssps_data_ssp5_2020 * adjfact_pop_ssp5_2020
shape$pop_ssps_data_ssp5_2030 <- shape$pop_ssps_data_ssp5_2030 * adjfact_pop_ssp5_2030
shape$pop_ssps_data_ssp5_2040 <- shape$pop_ssps_data_ssp5_2040 * adjfact_pop_ssp5_2040
shape$pop_ssps_data_ssp5_2050 <- shape$pop_ssps_data_ssp5_2050 * adjfact_pop_ssp5_2050

###

keep(shape, sure=T)

save.image("supporting_data/data_for_global_spline_v2.Rds")



## This R-script:
##      1) creates the gridded dataset of expenditure, CDDs/HDDs, and age, gender and sex which will be used to produce global gridded projection of future AC uptake and utilisation


rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

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

####

setwd(wd)

###


# gdp downscaled (SSPS)
gdp_ssps <- list.files(path="supporting_data/gdp_downscaled_ssps", recursive = T, pattern="tif", full.names = T)

gdp_ssps_data <- lapply(gdp_ssps, raster)
gdp_ssps_data <- split(gdp_ssps_data, rep(1:5, each=26))
gdp_ssps_data <- lapply(gdp_ssps_data, stack)
names(gdp_ssps_data[[1]]) <- paste0("GDP_SSP1_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[2]]) <- paste0("GDP_SSP2_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[3]]) <- paste0("GDP_SSP3_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[4]]) <- paste0("GDP_SSP4_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[5]]) <- paste0("GDP_SSP5_", seq(1850, 2100, by=10))

for (i in 1:5){
  gdp_ssps_data[[i]] <- raster::subset(gdp_ssps_data[[i]], 17:21)
}

# pop downscaled (SSPS)
pop_ssps <- list.files(path="supporting_data/pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- lapply(pop_ssps, stack)

for (i in 1:5){
  pop_ssps_data[[i]] <- raster::subset(pop_ssps_data[[i]], c(5+c(10*c(0:4))))
  
}

names(pop_ssps_data[[1]]) <- paste0("POP_SSP1_", seq(2010, 2050, by=10))
names(pop_ssps_data[[2]]) <- paste0("POP_SSP2_", seq(2010, 2050, by=10))
names(pop_ssps_data[[3]]) <- paste0("POP_SSP3_", seq(2010, 2050, by=10))
names(pop_ssps_data[[4]]) <- paste0("POP_SSP4_", seq(2010, 2050, by=10))
names(pop_ssps_data[[5]]) <- paste0("POP_SSP5_", seq(2010, 2050, by=10))

pop_ssps_data <- stack(brick(pop_ssps_data))

############

#gdp_ssps_data <- projectRaster(gdp_ssps_data, pop_ssps_data)

################

shape <- st_as_stars(pop_ssps_data[[1]])
shape <- st_as_sf(shape, na.rm=F)
st_crs(shape) <- 4326
shape <- st_transform(shape, 4326)

#add pop

pop_ssps_data_r <- exact_extract(pop_ssps_data, shape, "sum")

shape$POP_SSP2_2010 <- NULL

shape <- bind_cols(shape, pop_ssps_data_r)

# add macro-region

data("wrld_simpl")
wrld_simpl <- st_as_sf(wrld_simpl)

wrld_simpl$ISO3_n <- as.numeric(as.factor(wrld_simpl$ISO3))

raster_of_iso <- fasterize::fasterize(wrld_simpl, pop_ssps_data[[1]], field="ISO3_n", "first")

shape$ISO3_n <- exact_extract(raster_of_iso, shape, 'majority')

shape$id <- 1:nrow(shape)

shape_bk <- shape
wrld_simpl$geometry <- NULL
shape <- merge(shape, wrld_simpl, "ISO3_n", all.x=T)
shape$region <- countrycode(shape$ISO3, 'iso3c', 'region')
shape <- st_as_sf(shape)
shape <- shape[order(shape$id),]

#

adm1 <- read_sf("supporting_data/gadm_410-gpkg/gadm_410-levels.gpkg", layer = "ADM_1")
st_crs(adm1) <- 4326

rename_geometry <- function(g, name){
  current = attr(g, "sf_column")
  names(g)[names(g) == current] = name
  st_geometry(g) = name
  g
}

adm1 = rename_geometry(adm1, "geometry")

#adm1 <- st_make_valid(adm1)

###

gdp_ssps_data <- stack(gdp_ssps_data)
crs(gdp_ssps_data) <- crs(pop_ssps_data)

#

gdp_ssps_data_ex <- exact_extract(gdp_ssps_data, adm1, "sum", max_cells_in_memory = 93355200 )
pop_ssps_data_ex <- exact_extract(pop_ssps_data, adm1, "sum", max_cells_in_memory  = 93355200 )

#pop_ssps_data_ex <- pop_ssps_data_ex %>% mutate_all(~ifelse(.<100 | is.na(.), NA, .))

gdp_ssps_data_capita <- gdp_ssps_data_ex / pop_ssps_data_ex

gdp_ssps_data_capita = bind_cols(gdp_ssps_data_capita, adm1$GID_0)

# capper <- function(X){ifelse(X>quantile(X, 0.975, na.rm=T), quantile(X, 0.975, na.rm=T), ifelse(X<quantile(X, 0.025, na.rm=T), quantile(X, 0.025, na.rm=T), X))}
# 
# gdp_ssps_data_capita = gdp_ssps_data_capita %>% group_by(`...26`) %>% dplyr::mutate_all(., capper) 

gdp_ssps_data_capita <- gdp_ssps_data_capita %>% mutate_all(~ifelse(is.infinite(.) | .==0, NA, .))

gdp_ssps_data_capita <- gdp_ssps_data_capita %>% group_by(`...26`) %>%  mutate_all(~ifelse(is.na(.), median(., na.rm=T), .))

#tapply(gdp_ssps_data_capita$sum.GDP_SSP2_2050, gdp_ssps_data_capita$...26, summary)

gdp_ssps_data_capita$...26 = NULL

###

# difference between disposable income and consumption
#https://data.worldbank.org/indicator/NY.GNS.ICTR.ZS?view=map&year=2020

# library(wbstats)
# 
# savings_rate <- wb_data("NY.GNS.ICTR.ZS", mrnev=1)
# 
# gdp_ssps_data_capita <- bind_cols(gdp_ssps_data_capita, adm1$GID_0)
# 
# gdp_ssps_data_capita <- merge(gdp_ssps_data_capita, savings_rate, by.x="...26", by.y="iso3c", all.x=T)
# 
# gdp_ssps_data_capita$region <- countrycode::countrycode(gdp_ssps_data_capita$...26, 'iso3c', 'region')
# 
# gdp_ssps_data_capita <- gdp_ssps_data_capita %>% group_by(region) %>% mutate(NY.GNS.ICTR.ZS = ifelse(is.na(NY.GNS.ICTR.ZS), mean(NY.GNS.ICTR.ZS, na.rm=T), NY.GNS.ICTR.ZS))
# 
# gdp_ssps_data_capita$...26 <- NULL
# 
# gdp_ssps_data_capita <- gdp_ssps_data_capita[,c(1:10, 14)]
# 
# gdp_ssps_data_capita[,c(1:10)] <- gdp_ssps_data_capita[,c(1:10)] * (1 - (gdp_ssps_data_capita$NY.GNS.ICTR.ZS/100))
# 
# gdp_ssps_data_capita$NY.GNS.ICTR.ZS <- NULL

###############

gdp_ssps_data_capita <- st_as_sf(bind_cols(gdp_ssps_data_capita, adm1$geometry))

out <- list()

for (X in colnames(gdp_ssps_data_capita)[1:25]){
  out[[X]] <- fasterize(gdp_ssps_data_capita, raster_of_iso, X)
}

out <- stack(out)

ln_gdp_capita_sf <- exact_extract(out, shape, "median")

ln_gdp_capita_sf <- log(ln_gdp_capita_sf + 1)

ln_gdp_capita_sf <- as.data.frame(ln_gdp_capita_sf)
ln_gdp_capita_sf$id <- shape$id

###

shape_all <- merge(shape, ln_gdp_capita_sf, by="id")
shape_all <- st_as_sf(shape_all)

##################################

# Set directories
gldas <- "supporting_data/"

#cross calibrate CDDs/HDDs with GLDAS historical

cmip6_hist_cdd18_global  <-  c(paste0(gldas, 'gldas_0p25_deg_cdd_base_T_18C_1970_2019_ann.nc4', sep='')) 

cmip6_hist_cdd18_global  <- brick(cmip6_hist_cdd18_global, lvar=3, values=TRUE, level=1, 
                                  varname="cooling_degree_days_annual")

values(cmip6_hist_cdd18_global) <- ifelse(is.na(values(cmip6_hist_cdd18_global)), 0, values(cmip6_hist_cdd18_global))

CDD_glas_hist_sf <- exact_extract(cmip6_hist_cdd18_global, shape, 
                                  'mean', max_cells_in_memory = 93355200 )  

CDD_glas_hist_sf$geometry <- NULL

######################################################################	

#             CMIP6 Projections of Gridded Climate Data              #

######################################################################

##########
filter_gcm <- na.omit(readxl::read_xlsx("supporting_data/hausfather_hot_model_cmip6.xlsx", skip = 2)[,22])
filter_gcm <- filter_gcm$`TCR likely + ECS likely + GDDP`
###########

cmip6 <- "supporting_data/CDDs_18deg/hist"

cmip6_hist_cdd18_global_l <- list.files(cmip6, pattern = ".nc", full.names = T)
cmip6_hist_cdd18_global_l <- cmip6_hist_cdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_hist_cdd18_global_l)]

CDD_hist_sf_l <- list()

for(cmip6model in 1:length(cmip6_hist_cdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_hist_cdd18_global  <- cmip6_hist_cdd18_global_l[cmip6model] 
  
  cmip6_hist_cdd18_global  <- brick(cmip6_hist_cdd18_global, lvar=3, values=TRUE, level=1, 
                                    varname="tas")
  
  values(cmip6_hist_cdd18_global) <- ifelse(is.na(values(cmip6_hist_cdd18_global)), 0, values(cmip6_hist_cdd18_global))
  
  CDD_hist_sf <- exact_extract(cmip6_hist_cdd18_global, shape, 
                               'mean', max_cells_in_memory = 93355200 )  
  
  CDD_hist_sf$geometry <- NULL
  
  CDD_hist_sf_l[[cmip6model]] <- CDD_hist_sf
  
}

for(cmip6model in 1:length(cmip6_hist_cdd18_global_l)){
  
  model <- gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(CDD_hist_sf_l[[cmip6model]]) <- paste0(colnames(CDD_hist_sf_l[[cmip6model]]), "_", model)
  
}

CDD_hist_sf <- bind_cols(CDD_hist_sf_l)

CDD_hist_sf_mean <- sapply(split.default(CDD_hist_sf, ((1:ncol(CDD_hist_sf))-1)%/%20 + 1), rowMeans)
CDD_hist_sf_mean <- as.data.frame(CDD_hist_sf_mean)
colnames(CDD_hist_sf_mean) <- gsub(".nc", "", gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l)))


## CMIP6 Projections CDD - SSP2 RCP4.5
# Load netcdf file

cmip6 <- "supporting_data/CDDs_18deg/fut"

cmip6_245_cdd18_global_l <- list.files(cmip6, pattern = "245", full.names = T)
cmip6_245_cdd18_global_l <- cmip6_245_cdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_245_cdd18_global_l)]

CDD_245_sf_l <- list()

for(cmip6model in 1:length(cmip6_245_cdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_245_cdd18_global  <- cmip6_245_cdd18_global_l[cmip6model] 
  cmip6_245_cdd18_global  <- brick(cmip6_245_cdd18_global, lvar=3, values=TRUE, level=1, 
                                   varname="tas")
  
  
  values(cmip6_245_cdd18_global) <- ifelse(is.na(values(cmip6_245_cdd18_global)), 0, values(cmip6_245_cdd18_global))
  
  
  CDD_245_sf <- exact_extract(cmip6_245_cdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  CDD_245_sf$geometry <- NULL
  
  CDD_245_sf_l[[cmip6model]] <- CDD_245_sf
  
}

###########

for(cmip6model in 1:length(cmip6_245_cdd18_global_l)){
  
  model <- gsub("cdd_18_global_mid_245_", "", basename(cmip6_245_cdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(CDD_245_sf_l[[cmip6model]]) <- paste0(colnames(CDD_245_sf_l[[cmip6model]]), "_", model)
  
}

CDD_245_sf <- bind_cols(CDD_245_sf_l)

CDD_245_sf_mean <- sapply(split.default(CDD_245_sf, ((1:ncol(CDD_245_sf))-1)%/%36 + 1), rowMeans)
CDD_245_sf_mean <- as.data.frame(CDD_245_sf_mean)
colnames(CDD_245_sf_mean) <- gsub(".nc", "", gsub("cdd_18_global_mid_245_", "", basename(cmip6_245_cdd18_global_l)))

####

model <- gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l))
model <- gsub(".nc", "", model)

CDD_245_sf <- lapply(model, function(X){(CDD_245_sf %>% dplyr::select(ends_with(X))) - pull(CDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
CDD_245_sf <- bind_cols(CDD_245_sf)

CDD_245_sf <- rowMeans(CDD_glas_hist_sf[,c(31:50)], na.rm=T) + CDD_245_sf
CDD_245_sf[CDD_245_sf<0] <- 0

## CMIP6 Projections CDD - SSP5 RCP8.5
# Load netcdf file
cmip6 <- "supporting_data/CDDs_18deg/fut"

cmip6_585_cdd18_global_l <- list.files(cmip6, pattern = "585", full.names = T)
cmip6_585_cdd18_global_l <- cmip6_585_cdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_585_cdd18_global_l)]

CDD_585_sf_l <- list()

for(cmip6model in 1:length(cmip6_585_cdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_585_cdd18_global  <- cmip6_585_cdd18_global_l[cmip6model] 
  cmip6_585_cdd18_global  <- brick(cmip6_585_cdd18_global, lvar=3, values=TRUE, level=1, 
                                   varname="tas")
  
  
  values(cmip6_585_cdd18_global) <- ifelse(is.na(values(cmip6_585_cdd18_global)), 0, values(cmip6_585_cdd18_global))
  
  
  CDD_585_sf <- exact_extract(cmip6_585_cdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  CDD_585_sf$geometry <- NULL
  
  CDD_585_sf_l[[cmip6model]] <- CDD_585_sf
  
}

###########

for(cmip6model in 1:length(cmip6_585_cdd18_global_l)){
  
  model <- gsub("cdd_18_global_mid_585_", "", basename(cmip6_585_cdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(CDD_585_sf_l[[cmip6model]]) <- paste0(colnames(CDD_585_sf_l[[cmip6model]]), "_", model)
  
}

CDD_585_sf <- bind_cols(CDD_585_sf_l)

CDD_585_sf_mean <- sapply(split.default(CDD_585_sf, ((1:ncol(CDD_585_sf))-1)%/%36 + 1), rowMeans)
CDD_585_sf_mean <- as.data.frame(CDD_585_sf_mean)
colnames(CDD_585_sf_mean) <- gsub(".nc", "", gsub("cdd_18_global_mid_585_", "", basename(cmip6_585_cdd18_global_l)))

####

model <- gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l))
model <- gsub(".nc", "", model)

CDD_585_sf <- lapply(model, function(X){(CDD_585_sf %>% dplyr::select(ends_with(X))) - pull(CDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
CDD_585_sf <- bind_cols(CDD_585_sf)

CDD_585_sf <- rowMeans(CDD_glas_hist_sf[,c(31:50)], na.rm=T) + CDD_585_sf
CDD_585_sf[CDD_585_sf<0] <- 0

## CMIP6 Projections CDD - SSP1 RCP2.6
# Load netcdf file

cmip6 <- "supporting_data/CDDs_18deg/fut"

cmip6_126_cdd18_global_l <- list.files(cmip6, pattern = "126", full.names = T)
cmip6_126_cdd18_global_l <- cmip6_126_cdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_126_cdd18_global_l)]

CDD_126_sf_l <- list()

library(terra)

for(cmip6model in 1:length(cmip6_126_cdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_126_cdd18_global  <- cmip6_126_cdd18_global_l[cmip6model] 
  cmip6_126_cdd18_global  <- brick(rast(cmip6_126_cdd18_global))
  names(cmip6_126_cdd18_global) <- names(cmip6_245_cdd18_global)
  
  values(cmip6_126_cdd18_global) <- ifelse(is.na(values(cmip6_126_cdd18_global)), 0, values(cmip6_126_cdd18_global))
  
  
  CDD_126_sf <- exact_extract(cmip6_126_cdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  CDD_126_sf$geometry <- NULL
  
  CDD_126_sf_l[[cmip6model]] <- CDD_126_sf
  
}

###########

for(cmip6model in 1:length(cmip6_126_cdd18_global_l)){
  
  model <- gsub("cdd_18_global_mid_126_", "", basename(cmip6_126_cdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(CDD_126_sf_l[[cmip6model]]) <- paste0(colnames(CDD_126_sf_l[[cmip6model]]), "_", model)
  
}

CDD_126_sf <- bind_cols(CDD_126_sf_l)

CDD_126_sf_mean <- sapply(split.default(CDD_126_sf, ((1:ncol(CDD_126_sf))-1)%/%36 + 1), rowMeans)
CDD_126_sf_mean <- as.data.frame(CDD_126_sf_mean)
colnames(CDD_126_sf_mean) <- gsub(".nc", "", gsub("cdd_18_global_mid_126_", "", basename(cmip6_126_cdd18_global_l)))

####

model1 <- gsub("cdd_18_global_mid_126_", "", basename(cmip6_126_cdd18_global_l))
model1 <- gsub(".nc", "", model1)

model <- gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l))
model <- gsub(".nc", "", model)

model <- intersect(model, model1)

CDD_126_sf <- lapply(model, function(X){(CDD_126_sf %>% dplyr::select(ends_with(X))) - pull(CDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
CDD_126_sf <- bind_cols(CDD_126_sf)

CDD_126_sf <- rowMeans(CDD_glas_hist_sf[,c(31:50)], na.rm=T) + CDD_126_sf
CDD_126_sf[CDD_126_sf<0] <- 0

## CMIP6 Projections CDD - SSP1 RCP2.6
# Load netcdf file

cmip6 <- "supporting_data/CDDs_18deg/fut"

cmip6_370_cdd18_global_l <- list.files(cmip6, pattern = "370", full.names = T)
cmip6_370_cdd18_global_l <- cmip6_370_cdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_370_cdd18_global_l)]

CDD_370_sf_l <- list()

for(cmip6model in 1:length(cmip6_370_cdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_370_cdd18_global  <- cmip6_370_cdd18_global_l[cmip6model] 
  cmip6_370_cdd18_global  <- brick(rast(cmip6_370_cdd18_global))
  names(cmip6_370_cdd18_global) <- names(cmip6_245_cdd18_global)
  
  
  values(cmip6_370_cdd18_global) <- ifelse(is.na(values(cmip6_370_cdd18_global)), 0, values(cmip6_370_cdd18_global))
  
  
  CDD_370_sf <- exact_extract(cmip6_370_cdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  CDD_370_sf$geometry <- NULL
  
  CDD_370_sf_l[[cmip6model]] <- CDD_370_sf
  
}

###########

for(cmip6model in 1:length(cmip6_370_cdd18_global_l)){
  
  model <- gsub("cdd_18_global_mid_370_", "", basename(cmip6_370_cdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(CDD_370_sf_l[[cmip6model]]) <- paste0(colnames(CDD_370_sf_l[[cmip6model]]), "_", model)
  
}

CDD_370_sf <- bind_cols(CDD_370_sf_l)

CDD_370_sf_mean <- sapply(split.default(CDD_370_sf, ((1:ncol(CDD_370_sf))-1)%/%36 + 1), rowMeans)
CDD_370_sf_mean <- as.data.frame(CDD_370_sf_mean)
colnames(CDD_370_sf_mean) <- gsub(".nc", "", gsub("cdd_18_global_mid_370_", "", basename(cmip6_370_cdd18_global_l)))

####

model1 <- gsub("cdd_18_global_mid_370_", "", basename(cmip6_370_cdd18_global_l))
model1 <- gsub(".nc", "", model1)

model <- gsub("cdd_18_global_hist_", "", basename(cmip6_hist_cdd18_global_l))
model <- gsub(".nc", "", model)

model <- intersect(model, model1)

CDD_370_sf <- lapply(model, function(X){(CDD_370_sf %>% dplyr::select(ends_with(X))) - pull(CDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
CDD_370_sf <- bind_cols(CDD_370_sf)

CDD_370_sf <- rowMeans(CDD_glas_hist_sf[,c(31:50)], na.rm=T) + CDD_370_sf
CDD_370_sf[CDD_370_sf<0] <- 0

##############################

colnames(CDD_hist_sf) <- gsub("mean.X", "CDD_hist_", colnames(CDD_hist_sf))
colnames(CDD_245_sf) <- gsub("mean.X", "CDD_245_", colnames(CDD_245_sf))
colnames(CDD_585_sf) <- gsub("mean.X", "CDD_585_", colnames(CDD_585_sf))
colnames(CDD_126_sf) <- gsub("mean.X", "CDD_126_", colnames(CDD_126_sf))
colnames(CDD_370_sf) <- gsub("mean.X", "CDD_370_", colnames(CDD_370_sf))

##########

gldas <- "supporting_data/"

#cross calibrate HDDs/HDDs with GLDAS historical

cmip6_hist_hdd18_global  <-  c(paste0(gldas, 'gldas_0p25_deg_hdd_base_T_18C_1970_2019_ann.nc4', sep='')) 

cmip6_hist_hdd18_global  <- brick(cmip6_hist_hdd18_global, lvar=3, values=TRUE, level=1, 
                                  varname="heating_degree_days_annual")

values(cmip6_hist_hdd18_global) <- ifelse(is.na(values(cmip6_hist_hdd18_global)), 0, values(cmip6_hist_hdd18_global))

HDD_glas_hist_sf <- exact_extract(cmip6_hist_hdd18_global, shape, 
                                  'mean', max_cells_in_memory = 93355200 )  

HDD_glas_hist_sf$geometry <- NULL

###

######################################################################

#             CMIP6 Projections of Gridded Climate Data              #

######################################################################

cmip6 <- "supporting_data/HDDs_18deg/hist"

cmip6_hist_hdd18_global_l <- list.files(cmip6, pattern = ".nc", full.names = T)
cmip6_hist_hdd18_global_l <- cmip6_hist_hdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_hist_hdd18_global_l)]

HDD_hist_sf_l <- list()

for(cmip6model in 1:length(cmip6_hist_hdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_hist_hdd18_global  <- cmip6_hist_hdd18_global_l[cmip6model] 
  
  cmip6_hist_hdd18_global  <- brick(cmip6_hist_hdd18_global, lvar=3, values=TRUE, level=1, 
                                    varname="tas")
  
  values(cmip6_hist_hdd18_global) <- ifelse(is.na(values(cmip6_hist_hdd18_global)), 0, values(cmip6_hist_hdd18_global))
  
  HDD_hist_sf <- exact_extract(cmip6_hist_hdd18_global, shape, 
                               'mean', max_cells_in_memory = 93355200 )  
  
  HDD_hist_sf$geometry <- NULL
  
  HDD_hist_sf_l[[cmip6model]] <- HDD_hist_sf
  
}

for(cmip6model in 1:length(cmip6_hist_hdd18_global_l)){
  
  model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(HDD_hist_sf_l[[cmip6model]]) <- paste0(colnames(HDD_hist_sf_l[[cmip6model]]), "_", model)
  
}

HDD_hist_sf <- bind_cols(HDD_hist_sf_l)

HDD_hist_sf_mean <- sapply(split.default(HDD_hist_sf, ((1:ncol(HDD_hist_sf))-1)%/%20 + 1), rowMeans)
HDD_hist_sf_mean <- as.data.frame(HDD_hist_sf_mean)
colnames(HDD_hist_sf_mean) <- gsub(".nc", "", gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l)))


## CMIP6 Projections HDD - SSP2 RCP4.5
# Load netcdf file

cmip6 <- "supporting_data/HDDs_18deg/fut"

cmip6_245_hdd18_global_l <- list.files(cmip6, pattern = "245", full.names = T)
cmip6_245_hdd18_global_l <- cmip6_245_hdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_245_hdd18_global_l)]

HDD_245_sf_l <- list()

for(cmip6model in 1:length(cmip6_245_hdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_245_hdd18_global  <- cmip6_245_hdd18_global_l[cmip6model] 
  cmip6_245_hdd18_global  <- brick(cmip6_245_hdd18_global, lvar=3, values=TRUE, level=1, 
                                   varname="tas")
  
  
  values(cmip6_245_hdd18_global) <- ifelse(is.na(values(cmip6_245_hdd18_global)), 0, values(cmip6_245_hdd18_global))
  
  
  HDD_245_sf <- exact_extract(cmip6_245_hdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  HDD_245_sf$geometry <- NULL
  
  HDD_245_sf_l[[cmip6model]] <- HDD_245_sf
  
}

###########

for(cmip6model in 1:length(cmip6_245_hdd18_global_l)){
  
  model <- gsub("hdd_18_global_mid_245_", "", basename(cmip6_245_hdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(HDD_245_sf_l[[cmip6model]]) <- paste0(colnames(HDD_245_sf_l[[cmip6model]]), "_", model)
  
}

HDD_245_sf <- bind_cols(HDD_245_sf_l)

HDD_245_sf_mean <- sapply(split.default(HDD_245_sf, ((1:ncol(HDD_245_sf))-1)%/%36 + 1), rowMeans)
HDD_245_sf_mean <- as.data.frame(HDD_245_sf_mean)
colnames(HDD_245_sf_mean) <- gsub(".nc", "", gsub("hdd_18_global_mid_245_", "", basename(cmip6_245_hdd18_global_l)))

####

model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l))
model <- gsub(".nc", "", model)

HDD_245_sf <- lapply(model, function(X){(HDD_245_sf %>% dplyr::select(ends_with(X))) - pull(HDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
HDD_245_sf <- bind_cols(HDD_245_sf)

HDD_245_sf <- rowMeans(HDD_glas_hist_sf[,c(31:50)], na.rm=T) + HDD_245_sf
HDD_245_sf[HDD_245_sf<0] <- 0

## CMIP6 Projections HDD - SSP5 RCP8.5
# Load netcdf file
cmip6 <- "supporting_data/HDDs_18deg/fut"

cmip6_585_hdd18_global_l <- list.files(cmip6, pattern = "585", full.names = T)
cmip6_585_hdd18_global_l <- cmip6_585_hdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_585_hdd18_global_l)]

HDD_585_sf_l <- list()

for(cmip6model in 1:length(cmip6_585_hdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_585_hdd18_global  <- cmip6_585_hdd18_global_l[cmip6model] 
  cmip6_585_hdd18_global  <- brick(cmip6_585_hdd18_global, lvar=3, values=TRUE, level=1, 
                                   varname="tas")
  
  
  values(cmip6_585_hdd18_global) <- ifelse(is.na(values(cmip6_585_hdd18_global)), 0, values(cmip6_585_hdd18_global))
  
  
  HDD_585_sf <- exact_extract(cmip6_585_hdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  HDD_585_sf$geometry <- NULL
  
  HDD_585_sf_l[[cmip6model]] <- HDD_585_sf
  
}

###########

for(cmip6model in 1:length(cmip6_585_hdd18_global_l)){
  
  model <- gsub("hdd_18_global_mid_585_", "", basename(cmip6_585_hdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(HDD_585_sf_l[[cmip6model]]) <- paste0(colnames(HDD_585_sf_l[[cmip6model]]), "_", model)
  
}

HDD_585_sf <- bind_cols(HDD_585_sf_l)

HDD_585_sf_mean <- sapply(split.default(HDD_585_sf, ((1:ncol(HDD_585_sf))-1)%/%36 + 1), rowMeans)
HDD_585_sf_mean <- as.data.frame(HDD_585_sf_mean)
colnames(HDD_585_sf_mean) <- gsub(".nc", "", gsub("hdd_18_global_mid_585_", "", basename(cmip6_585_hdd18_global_l)))

####

model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l))
model <- gsub(".nc", "", model)

HDD_585_sf <- lapply(model, function(X){(HDD_585_sf %>% dplyr::select(ends_with(X))) - pull(HDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
HDD_585_sf <- bind_cols(HDD_585_sf)

HDD_585_sf <- rowMeans(HDD_glas_hist_sf[,c(31:50)], na.rm=T) + HDD_585_sf
HDD_585_sf[HDD_585_sf<0] <- 0

## CMIP6 Projections HDD - SSP1 RCP2.6
# Load netcdf file

cmip6 <- "supporting_data/HDDs_18deg/fut"

cmip6_126_hdd18_global_l <- list.files(cmip6, pattern = "126", full.names = T)
cmip6_126_hdd18_global_l <- cmip6_126_hdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_126_hdd18_global_l)]

HDD_126_sf_l <- list()

for(cmip6model in 1:length(cmip6_126_hdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_126_hdd18_global  <- cmip6_126_hdd18_global_l[cmip6model] 
  cmip6_126_hdd18_global  <- brick(rast(cmip6_126_hdd18_global))
  names(cmip6_126_hdd18_global) <- names(cmip6_245_hdd18_global)
  
  values(cmip6_126_hdd18_global) <- ifelse(is.na(values(cmip6_126_hdd18_global)), 0, values(cmip6_126_hdd18_global))
  
  
  HDD_126_sf <- exact_extract(cmip6_126_hdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  HDD_126_sf$geometry <- NULL
  
  HDD_126_sf_l[[cmip6model]] <- HDD_126_sf
  
}

###########

for(cmip6model in 1:length(cmip6_126_hdd18_global_l)){
  
  model <- gsub("hdd_18_global_mid_126_", "", basename(cmip6_126_hdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(HDD_126_sf_l[[cmip6model]]) <- paste0(colnames(HDD_126_sf_l[[cmip6model]]), "_", model)
  
}

HDD_126_sf <- bind_cols(HDD_126_sf_l)

HDD_126_sf_mean <- sapply(split.default(HDD_126_sf, ((1:ncol(HDD_126_sf))-1)%/%36 + 1), rowMeans)
HDD_126_sf_mean <- as.data.frame(HDD_126_sf_mean)
colnames(HDD_126_sf_mean) <- gsub(".nc", "", gsub("hdd_18_global_mid_126_", "", basename(cmip6_126_hdd18_global_l)))

####

model1 <- gsub("hdd_18_global_mid_126_", "", basename(cmip6_126_hdd18_global_l))
model1 <- gsub(".nc", "", model1)

model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l))
model <- gsub(".nc", "", model)

model <- intersect(model, model1)

HDD_126_sf <- lapply(model, function(X){(HDD_126_sf %>% dplyr::select(ends_with(X))) - pull(HDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
HDD_126_sf <- bind_cols(HDD_126_sf)

HDD_126_sf <- rowMeans(HDD_glas_hist_sf[,c(31:50)], na.rm=T) + HDD_126_sf
HDD_126_sf[HDD_126_sf<0] <- 0

## CMIP6 Projections HDD - SSP1 RCP2.6
# Load netcdf file

cmip6 <- "supporting_data/HDDs_18deg/fut"

cmip6_370_hdd18_global_l <- list.files(cmip6, pattern = "370", full.names = T)
cmip6_370_hdd18_global_l <- cmip6_370_hdd18_global_l[grepl(paste(filter_gcm, collapse = "|"), cmip6_370_hdd18_global_l)]

HDD_370_sf_l <- list()

for(cmip6model in 1:length(cmip6_370_hdd18_global_l)){
  
  print(cmip6model)
  
  cmip6_370_hdd18_global  <- cmip6_370_hdd18_global_l[cmip6model] 
  cmip6_370_hdd18_global  <- brick(rast(cmip6_370_hdd18_global))
  names(cmip6_370_hdd18_global) <- names(cmip6_245_hdd18_global)
  
  values(cmip6_370_hdd18_global) <- ifelse(is.na(values(cmip6_370_hdd18_global)), 0, values(cmip6_370_hdd18_global))
  
  
  HDD_370_sf <- exact_extract(cmip6_370_hdd18_global, shape, 
                              'mean', max_cells_in_memory = 93355200 )   #  max_cells_in_memory = 3e+07 is default
  
  HDD_370_sf$geometry <- NULL
  
  HDD_370_sf_l[[cmip6model]] <- HDD_370_sf
  
}

###########

for(cmip6model in 1:length(cmip6_370_hdd18_global_l)){
  
  model <- gsub("hdd_18_global_mid_370_", "", basename(cmip6_370_hdd18_global_l)[cmip6model])
  model <- gsub(".nc", "", model)
  
  colnames(HDD_370_sf_l[[cmip6model]]) <- paste0(colnames(HDD_370_sf_l[[cmip6model]]), "_", model)
  
}

HDD_370_sf <- bind_cols(HDD_370_sf_l)

HDD_370_sf_mean <- sapply(split.default(HDD_370_sf, ((1:ncol(HDD_370_sf))-1)%/%36 + 1), rowMeans)
HDD_370_sf_mean <- as.data.frame(HDD_370_sf_mean)
colnames(HDD_370_sf_mean) <- gsub(".nc", "", gsub("hdd_18_global_mid_370_", "", basename(cmip6_370_hdd18_global_l)))

####

model1 <- gsub("hdd_18_global_mid_370_", "", basename(cmip6_370_hdd18_global_l))
model1 <- gsub(".nc", "", model1)

model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l))
model <- gsub(".nc", "", model)

model <- intersect(model, model1)

HDD_370_sf <- lapply(model, function(X){(HDD_370_sf %>% dplyr::select(ends_with(X))) - pull(HDD_hist_sf_mean %>% dplyr::select(ends_with(X)))})
HDD_370_sf <- bind_cols(HDD_370_sf)

HDD_370_sf <- rowMeans(HDD_glas_hist_sf[,c(31:50)], na.rm=T) + HDD_370_sf
HDD_370_sf[HDD_370_sf<0] <- 0

##############################

colnames(HDD_hist_sf) <- gsub("mean.X", "HDD_hist_", colnames(HDD_hist_sf))
colnames(HDD_245_sf) <- gsub("mean.X", "HDD_245_", colnames(HDD_245_sf))
colnames(HDD_585_sf) <- gsub("mean.X", "HDD_585_", colnames(HDD_585_sf))
colnames(HDD_126_sf) <- gsub("mean.X", "HDD_126_", colnames(HDD_126_sf))
colnames(HDD_370_sf) <- gsub("mean.X", "HDD_370_", colnames(HDD_370_sf))

###########
###########

CDD_126_sf$id <- shape$id
CDD_370_sf$id <- shape$id
CDD_245_sf$id <- shape$id
CDD_585_sf$id <- shape$id
HDD_126_sf$id <- shape$id
HDD_370_sf$id <- shape$id
HDD_245_sf$id <- shape$id
HDD_585_sf$id <- shape$id

###########

shape_all_bind <- Reduce(function(x,y) merge(x,y,by="id",all=TRUE) ,list(shape_all, CDD_126_sf, CDD_370_sf, CDD_245_sf, CDD_585_sf, HDD_126_sf, HDD_370_sf, HDD_245_sf, HDD_585_sf))

shape_all_bind <- dplyr::mutate_at(shape_all_bind, vars(contains('CDD') | contains('HDD')), ~ifelse(.<0, 0, .))

shape_all <- shape_all_bind

##################################

# urbanisation

# urban share downscaled (SSPS)
urban_ssps <- list.files(path="urban_share_downscaled_ssps/UrbanFraction_1_8_dgr_NETCDF_Projections_SSPs1-5_2010-2100_v1/", recursive = T, pattern="nc", full.names = T)

urban_ssps_data <- lapply(urban_ssps, raster)
urban_ssps_data <- split(urban_ssps_data, rep(1:5, each=10))
urban_ssps_data <- lapply(urban_ssps_data, stack)

names(urban_ssps_data[[1]]) <- paste0("URB_SSP1_", seq(2010, 2100, by=10))
names(urban_ssps_data[[2]]) <- paste0("URB_SSP2_", seq(2010, 2100, by=10))
names(urban_ssps_data[[3]]) <- paste0("URB_SSP3_", seq(2010, 2100, by=10))
names(urban_ssps_data[[4]]) <- paste0("URB_SSP4_", seq(2010, 2100, by=10))
names(urban_ssps_data[[5]]) <- paste0("URB_SSP5_", seq(2010, 2100, by=10))

#

urban_ssps_data_out <- list()
pop_ssps_data_w <- pop_ssps_data
for (i in 1:length(urban_ssps_data)){
  urban_ssps_data_out[[i]] <- exact_extract(urban_ssps_data[[i]], shape, fun="weighted_mean", weights=pop_ssps_data_w[[i]], max_cells_in_memory = 1e9)
}

urban_ssps_data_out_c <- bind_cols(urban_ssps_data_out)

colnames(urban_ssps_data_out_c) <- c(paste0("URB_SSP1_", seq(2010, 2100, by=10)), paste0("URB_SSP2_", seq(2010, 2100, by=10)), paste0("URB_SSP3_", seq(2010, 2100, by=10)), paste0("URB_SSP4_", seq(2010, 2100, by=10)), paste0("URB_SSP5_", seq(2010, 2100, by=10)))

urban_ssps_data_out_c$id <- shape$id

shape_all <- merge(shape_all, urban_ssps_data_out_c, "id")

##############################

# add age gender education

# Gender head (national SSPs)

# pop features (SSPS)
pop_features <- read.csv("pop_features_ssps/SspDb_country_data_2013-06-12.csv")

pop_features$SCENARIO <- substr(pop_features$SCENARIO, 1, 4)

pop_features <- pop_features %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

pop_features <- dplyr::select(pop_features, 1,2,3,4,18, 20, 22, 24, 26, 28, 30, 32, 34, 36)

pop_features_all <- filter(pop_features, VARIABLE=="Population") %>% ungroup()
pop_features_all <- dplyr::select(pop_features_all, c(2,3, 5:14))

pop_features_all <- reshape2::melt(pop_features_all, c(1,2))
pop_features_all$variable <- as.numeric(as.character(gsub("X", "", pop_features_all$variable )))

pop_features_all <- group_by(pop_features_all, SCENARIO, REGION, variable) %>% dplyr::summarise(value=first(value))

pop_features_all_before_melt <- pop_features_all

pop_features_all <- pivot_wider(
  pop_features_all,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "pop_",
  values_fn = first
)

pop_features_gender <- filter(pop_features, VARIABLE=="Population|Female")  %>% ungroup()
pop_features_gender <- dplyr::select(pop_features_gender, c(2,3, 5:14))

pop_features_gender <- reshape2::melt(pop_features_gender, c(1,2))
pop_features_gender$variable <- as.numeric(as.character(gsub("X", "", pop_features_gender$variable )))

pop_features_gender <- pop_features_gender[order( pop_features_gender$SCENARIO, pop_features_gender$REGION, pop_features_gender$variable ),]

pop_features_all_before_melt <- filter(pop_features_all_before_melt, REGION %in% pop_features_gender$REGION)

pop_features_gender$value <- pop_features_gender$value / pop_features_all_before_melt$value

pop_features_gender <- pivot_wider(
  pop_features_gender,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "women_"
)

shape_all <- merge(shape_all, pop_features_gender, by.x="ISO3", by.y="REGION", all.x=T)

# Age  (national SSPs)
pop_features_age <-filter(pop_features, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))
pop_features_age <- dplyr::select(pop_features_age, c(2:4, 5:14))

pop_features_age$VARIABLE <- gsub("Population|Female|Aged", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE <- gsub("Population|Male|Aged", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE  <- gsub("\\|", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE  <- gsub("\\+", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE <- sapply(strsplit(pop_features_age$VARIABLE, split = "-", fixed = TRUE), function(k) mean(as.numeric(k)))

pop_features_age <- dplyr::select(pop_features_age, c(1, 2, 3, 4:8))

colnames(pop_features_age)[4:8] <- paste0("X", seq(2010, 2050, 10))

pop_features_age <- pop_features_age %>% group_by(SCENARIO, REGION) %>% dplyr::summarise(X2010=sum((X2010/sum(X2010))*VARIABLE), X2020=sum((X2020/sum(X2020))*VARIABLE), X2030=sum((X2030/sum(X2030))*VARIABLE), X2040=sum((X2040/sum(X2040))*VARIABLE), X2050=sum((X2050/sum(X2050))*VARIABLE))

pop_features_age <- reshape2::melt(pop_features_age, c(1, 2))
pop_features_age$variable <- as.numeric(as.character(gsub("X", "", pop_features_age$variable )))

pop_features_age <- pivot_wider(
  pop_features_age,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "age_"
)

shape_all <- merge(shape_all, pop_features_age, by.x="ISO3", by.y="REGION", all.x=T)

# Edu (national SSPs)

pop_features_edu <-filter(pop_features, grepl("Education",VARIABLE))
pop_features_edu <- dplyr::select(pop_features_edu, c(2:4, 5:9))
pop_features_edu$VARIABLE <- gsub('.*\\|', "", pop_features_edu$VARIABLE)

pop_features_edu <- group_by(pop_features_edu, SCENARIO, VARIABLE, REGION) %>% dplyr::summarise_all(., "sum")

pop_features_edu <- group_by(pop_features_edu, SCENARIO, REGION) %>% dplyr::mutate(X2010=X2010/sum(X2010), X2020=X2020/sum(X2020), X2030=X2030/sum(X2030), X2040=X2040/sum(X2040), X2050=X2050/sum(X2050))

colnames(pop_features_edu)[4:8] <- paste0("X", seq(2010, 2050, 10))

pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="No Education"] <- "0"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Primary Education"] <- "1"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Secondary Education"] <-"2"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Tertiary Education"] <- "3"
pop_features_edu$VARIABLE <- as.numeric(pop_features_edu$VARIABLE)

pop_features_edu <- pop_features_edu %>% group_by(SCENARIO, REGION) %>% dplyr::summarise(X2010=sum((X2010/sum(X2010))*VARIABLE), X2020=sum((X2020/sum(X2020))*VARIABLE), X2030=sum((X2030/sum(X2030))*VARIABLE), X2040=sum((X2040/sum(X2040))*VARIABLE), X2050=sum((X2050/sum(X2050))*VARIABLE))

pop_features_edu <- reshape2::melt(pop_features_edu, c(1, 2))
pop_features_edu$variable <- as.numeric(as.character(gsub("X", "", pop_features_edu$variable )))

pop_features_edu <- pivot_wider(
  pop_features_edu,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "edu_"
)

shape_all <- merge(shape_all, pop_features_edu, by.x="ISO3", by.y="REGION", all.x=T)

##################################

shape <- st_as_sf(shape_all)

colnames(shape) <- gsub("median.sum.", "", colnames(shape))

colnames(shape)[4] <- "pop_ssps_data_hist"
v <- expand.grid("pop_ssps_data", "_", seq(2010, 2050, 10), "_", c("ssp1", "ssp2", "ssp3", "ssp4", "ssp5"))
colnames(shape)[5:29] <- paste0(v$Var1, v$Var2, v$Var5, v$Var4, v$Var3)

####

# add humidity

# load(paste0(wd, "/supporting_data/data_for_global_spline_v2.Rds"))

cmip6 <- "L:/falchetta/humidity_cmip"

cmip6_hist_hurs_global_l <- list.files(cmip6, pattern = ".nc", full.names = T, recursive = F)
cmip6_hist_hurs_global_l <- cmip6_hist_hurs_global_l[grepl("hist", cmip6_hist_hurs_global_l)]

cmip6_hist_hurs_global_l <- stack(lapply(cmip6_hist_hurs_global_l, raster, band=10))
cmip6_hist_hurs_global_l <- rotate(cmip6_hist_hurs_global_l)

cmip6_hist_hurs_global_l <- stackApply(cmip6_hist_hurs_global_l, 1, mean, na.rm = T)

shape <- st_as_sf(shape)

hurs_hist_sf <- exact_extract(cmip6_hist_hurs_global_l, shape, 
                              'mean', max_cells_in_memory = 93355200 )  

hurs_hist_sf <- as.data.frame(hurs_hist_sf)

colnames(hurs_hist_sf) <- c("hurs_2010")

namess <- do.call(paste0, expand.grid(colnames(hurs_hist_sf), "_", filter_gcm))
hurs_hist_sf <- bind_cols(rep(hurs_hist_sf, 14))
colnames(hurs_hist_sf) <- namess

###

cmip6 <- "L:/falchetta/humidity_cmip"

cmip6_126_hurs_global_l <- list.files(cmip6, pattern = ".nc", full.names = T, recursive = F)
cmip6_126_hurs_global_l <- cmip6_126_hurs_global_l[grepl("126", cmip6_126_hurs_global_l)]
cmip6_126_hurs_global_l <- cmip6_126_hurs_global_l[!grepl("json", cmip6_126_hurs_global_l)]
cmip6_126_hurs_global_l <- stack(lapply(cmip6_126_hurs_global_l, raster, band=10))
cmip6_126_hurs_global_l <- rotate(cmip6_126_hurs_global_l)
cmip6_126_hurs_global_l <- stackApply(cmip6_126_hurs_global_l, rep(2015:2055, each=12), mean, na.rm = T)
cmip6_126_hurs_global_l <- stackApply(cmip6_126_hurs_global_l, c(1, rep(1:4, each=10)), mean, na.rm = T)

hurs_126_sf <- exact_extract(cmip6_126_hurs_global_l, shape, 
                             'mean', max_cells_in_memory = 93355200 )  

colnames(hurs_126_sf) <- c("hurs_126_2020", "hurs_126_2030", "hurs_126_2040", "hurs_126_2050")

namess <- do.call(paste0, expand.grid(colnames(hurs_126_sf), "_", filter_gcm))
hurs_126_sf <- bind_cols(rep(hurs_126_sf, 14))
colnames(hurs_126_sf) <- namess

###

cmip6 <- "L:/falchetta/humidity_cmip"

cmip6_245_hurs_global_l <- list.files(cmip6, pattern = ".nc", full.names = T, recursive = F)
cmip6_245_hurs_global_l <- cmip6_245_hurs_global_l[grepl("245", cmip6_245_hurs_global_l)]
cmip6_245_hurs_global_l <- cmip6_245_hurs_global_l[!grepl("json", cmip6_245_hurs_global_l)]
cmip6_245_hurs_global_l <- stack(lapply(cmip6_245_hurs_global_l, raster, band=10))

cmip6_245_hurs_global_l <- rotate(cmip6_245_hurs_global_l)

cmip6_245_hurs_global_l <- stackApply(cmip6_245_hurs_global_l, rep(2015:2055, each=12), mean, na.rm = T)
cmip6_245_hurs_global_l <- stackApply(cmip6_245_hurs_global_l, c(1, rep(1:4, each=10)), mean, na.rm = T)

hurs_245_sf <- exact_extract(cmip6_245_hurs_global_l, shape, 
                             'mean', max_cells_in_memory = 93355200 )  

colnames(hurs_245_sf) <- c("hurs_245_2020", "hurs_245_2030", "hurs_245_2040", "hurs_245_2050")

namess <- do.call(paste0, expand.grid(colnames(hurs_245_sf), "_", filter_gcm))
hurs_245_sf <- bind_cols(rep(hurs_245_sf, 14))
colnames(hurs_245_sf) <- namess

###

cmip6 <- "L:/falchetta/humidity_cmip"

cmip6_370_hurs_global_l <- list.files(cmip6, pattern = ".nc", full.names = T, recursive = F)
cmip6_370_hurs_global_l <- cmip6_370_hurs_global_l[grepl("370", cmip6_370_hurs_global_l)]
cmip6_370_hurs_global_l <- cmip6_370_hurs_global_l[!grepl("json", cmip6_370_hurs_global_l)]
cmip6_370_hurs_global_l <- stack(lapply(cmip6_370_hurs_global_l, raster, band=10))
cmip6_370_hurs_global_l <- rotate(cmip6_370_hurs_global_l)
cmip6_370_hurs_global_l <- stackApply(cmip6_370_hurs_global_l, rep(2015:2055, each=12), mean, na.rm = T)
cmip6_370_hurs_global_l <- stackApply(cmip6_370_hurs_global_l, c(1, rep(1:4, each=10)), mean, na.rm = T)

hurs_370_sf <- exact_extract(cmip6_370_hurs_global_l, shape, 
                             'mean', max_cells_in_memory = 93355200 )  

colnames(hurs_370_sf) <- c("hurs_370_2020", "hurs_370_2030", "hurs_370_2040", "hurs_370_2050")

namess <- do.call(paste0, expand.grid(colnames(hurs_370_sf), "_", filter_gcm))
hurs_370_sf <- bind_cols(rep(hurs_370_sf, 14))
colnames(hurs_370_sf) <- namess


###

cmip6 <- "L:/falchetta/humidity_cmip"

cmip6_585_hurs_global_l <- list.files(cmip6, pattern = ".nc", full.names = T, recursive = F)
cmip6_585_hurs_global_l <- cmip6_585_hurs_global_l[grepl("585", cmip6_585_hurs_global_l)]
cmip6_585_hurs_global_l <- cmip6_585_hurs_global_l[!grepl("json", cmip6_585_hurs_global_l)]
cmip6_585_hurs_global_l <- stack(lapply(cmip6_585_hurs_global_l, raster, band=10))
cmip6_585_hurs_global_l <- rotate(cmip6_585_hurs_global_l)
cmip6_585_hurs_global_l <- stackApply(cmip6_585_hurs_global_l, rep(2015:2055, each=12), mean, na.rm = T)
cmip6_585_hurs_global_l <- stackApply(cmip6_585_hurs_global_l, c(1, rep(1:4, each=10)), mean, na.rm = T)

hurs_585_sf <- exact_extract(cmip6_585_hurs_global_l, shape, 
                             'mean', max_cells_in_memory = 93355200 )  

colnames(hurs_585_sf) <- c("hurs_585_2020", "hurs_585_2030", "hurs_585_2040", "hurs_585_2050")

namess <- do.call(paste0, expand.grid(colnames(hurs_585_sf), "_", filter_gcm))
hurs_585_sf <- bind_cols(rep(hurs_585_sf, 14))
colnames(hurs_585_sf) <- namess

hurs_585_sf

####

shape <- bind_cols(shape, hurs_hist_sf, hurs_126_sf, hurs_245_sf, hurs_370_sf, hurs_585_sf)

####

# add prices

setwd(wd)

library(readxl)
prices_ely <- read_xlsx("supporting_data/ely_prices.xlsx")

prices_ely$ctr <- prices_ely$`Country code`

# convert to 2011 PPP

library(wbstats)

ppp <- wb(indicator = "PA.NUS.PPP")
ppp <- group_by(ppp, iso2c) %>% dplyr::summarise(value=value[date==2021]/value[date==2011])

prices_ely <- merge(prices_ely, ppp, by.x="ctr", by.y="iso2c")

prices_ely$elyprc <- prices_ely$`Average price of 1KW/h (USD)` * prices_ely$value

prices_ely <- dplyr::select(prices_ely, ctr, elyprc)

shape <- merge(shape, prices_ely, by.x="ISO2", by.y="ctr", all.x=T)

###

# NA inputation

library(imputeTS)
library(missRanger)

c <- st_coordinates(st_centroid(shape))
bb <- bind_cols(shape$elyprc, as.data.frame(c))
bb <-missRanger::missRanger(bb, seed = 1234, data_only=T)

shape$elyprc <- bb$...1

shape$elyprc[is.na(shape$ISO3)] <- NA

###

shape_nogeo <- shape %>% st_set_geometry(NULL)

which_nas <- lapply(1:ncol(shape_nogeo), function(X){sum(is.na(shape_nogeo[!is.na(shape_nogeo$ISO3),X]))})
names(which_nas) <- colnames(shape_nogeo)

###

save(shape, file=paste0(wd, "/supporting_data/data_for_global_spline_v2.Rds"))

##########################
##########################
##########################

model1 <- gsub("hdd_18_global_mid_370_", "", basename(cmip6_370_hdd18_global_l))
model1 <- gsub(".nc", "", model1)

model <- gsub("hdd_18_global_hist_", "", basename(cmip6_hist_hdd18_global_l))
model <- gsub(".nc", "", model)

models_list <- intersect(model, model1)

save(models_list, file=paste0(wd, "/supporting_data/models_list.Rds"))

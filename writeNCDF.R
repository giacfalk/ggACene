rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

library(raster)
library(fasterize)
library(sf)

setwd(wd)

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

# write gridded data

shape_ac <- dplyr::select(shape_ac, starts_with("SSP"))
combs <- colnames(shape_ac)

r <- list.files(path="supporting_data/pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
r <- lapply(r, stack)
r <- r[[2]][[1]]

output_ac <- lapply(combs[1:20], function(X){fasterize(shape_ac, r, X)})
output_ac <- stack(output_ac)
names(output_ac) <- combs[1:20]


########################
########################

library(ncdf4)

filename <- "output_data/ncdf/ac_penetration_ssp1_ssp2_ssp3_ssp5_2010_2050.nc"

xvals <- unique(values(init(output_ac, "x")))
yvals <- unique(values(init(output_ac, "y")))
nx <- length(xvals)
ny <- length(yvals)
lon <- ncdim_def("longitude", "degrees_east", xvals)
lat <- ncdim_def("latitude", "degrees_north", yvals)

mv <- -999

time <- ncdim_def(name = "Time", 
                  units = "years", 
                  vals = c(2010, 2020, 2030, 2040, 2050), 
                  unlim = TRUE,
                  longname = "Year")

scen <- ncdim_def(name = "SSP", 
                  units = "scenario", 
                  vals = c(1, 2, 3, 5), 
                  unlim = TRUE,
                  longname = "SSP scenario")

var_prec <- ncvar_def(name = "ac_penetration",
                      units = "share_of_households",
                      dim = list(lon, lat, scen, time),
                      longname = "ac_penetration_rate",
                      missval = mv,
                      compression = 9)

ncout <- nc_create(filename, list(var_prec), force_v4 = TRUE)

ncatt_put(ncout, 0, "Title", "AC_penetration")
ncatt_put(ncout, 0, "Source", "Falchetta et al. (2023)")
ncatt_put(ncout, 0, "References", "Falchetta et al. (2023)")
ncatt_put(ncout, 0, "Created on", date())

for (i in c(1:5)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(output_ac[[i]]), 
            start = c(1, 1, 1, match(i, c(1:5))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(6:10)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(output_ac[[i]]), 
            start = c(1, 1, 2, match(i, c(6:10))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(11:15)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(output_ac[[i]]), 
            start = c(1, 1, 3, match(i, c(11:15))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(16:20)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(output_ac[[i]]), 
            start = c(1, 1, 4, match(i, c(16:20))), 
            count = c(-1, -1, 1, 1))
}

nc_close(ncout)

#####################################
#####################################

shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")
shape_ely_diff <- dplyr::select(shape_ely_diff, contains("total"))
combs <- colnames(shape_ely_diff)

r <- list.files(path="supporting_data/pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
r <- lapply(r, stack)
r <- r[[2]][[1]]

shape_ely_diff <- lapply(combs[1:20], function(X){fasterize(shape_ely_diff, r, X)})
shape_ely_diff <- stack(shape_ely_diff)/1e3
names(shape_ely_diff) <- combs[1:20]


####

filename <- "output_data/ncdf/ac_TWh_ssp1_ssp2_ssp3_ssp5_2010_2050.nc"

xvals <- unique(values(init(shape_ely_diff, "x")))
yvals <- unique(values(init(shape_ely_diff, "y")))
nx <- length(xvals)
ny <- length(yvals)
lon <- ncdim_def("longitude", "degrees_east", xvals)
lat <- ncdim_def("latitude", "degrees_north", yvals)

mv <- -999

time <- ncdim_def(name = "Time", 
                  units = "years", 
                  vals = c(2010, 2020, 2030, 2040, 2050), 
                  unlim = TRUE,
                  longname = "Year")

scen <- ncdim_def(name = "SSP", 
                  units = "scenario", 
                  vals = c(1, 2, 3, 5), 
                  unlim = TRUE,
                  longname = "SSP scenario")

var_prec <- ncvar_def(name = "ac_electricity_consumption",
                      units = "TWh/yr",
                      dim = list(lon, lat, scen, time),
                      longname = "ac_electricity_consumption",
                      missval = mv,
                      compression = 9)

ncout <- nc_create(filename, list(var_prec), force_v4 = TRUE)

ncatt_put(ncout, 0, "Title", "AC_electricity_consumption")
ncatt_put(ncout, 0, "Source", "Falchetta et al. (2023)")
ncatt_put(ncout, 0, "References", "Falchetta et al. (2023)")
ncatt_put(ncout, 0, "Created on", date())

for (i in c(1:5)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 1, match(i, c(1:5))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(6:10)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 2, match(i, c(6:10))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(11:15)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 3, match(i, c(11:15))), 
            count = c(-1, -1, 1, 1))
}


for (i in c(16:20)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 4, match(i, c(16:20))), 
            count = c(-1, -1, 1, 1))
}


nc_close(ncout)

######

shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")
shape_ely_diff <- dplyr::select(shape_ely_diff, contains("pop_ssps"))
combs <- colnames(shape_ely_diff)

r <- list.files(path="supporting_data/pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
r <- lapply(r, stack)
r <- r[[2]][[1]]

shape_ely_diff <- lapply(combs[1:17], function(X){fasterize(shape_ely_diff, r, X)})
shape_ely_diff <- stack(shape_ely_diff)
names(shape_ely_diff) <- combs[1:17]
shape_ely_diff <- stack(shape_ely_diff[[1:5]], shape_ely_diff[[1]], shape_ely_diff[[6:9]], shape_ely_diff[[1]], shape_ely_diff[[10:13]], shape_ely_diff[[1]], shape_ely_diff[[14:17]])
names(shape_ely_diff)[1] <- "pop_ssps_data_ssp1_2010"
names(shape_ely_diff)[6] <- "pop_ssps_data_ssp2_2010"
names(shape_ely_diff)[10] <- "pop_ssps_data_ssp3_2010"
names(shape_ely_diff)[14] <- "pop_ssps_data_ssp5_2010"

filename <- "output_data/ncdf/pop_ssp1_ssp2_ssp3_ssp5_2010_2050.nc"

xvals <- unique(values(init(shape_ely_diff, "x")))
yvals <- unique(values(init(shape_ely_diff, "y")))
nx <- length(xvals)
ny <- length(yvals)
lon <- ncdim_def("longitude", "degrees_east", xvals)
lat <- ncdim_def("latitude", "degrees_north", yvals)

mv <- -999

time <- ncdim_def(name = "Time", 
                  units = "years", 
                  vals = c(2010, 2020, 2030, 2040, 2050), 
                  unlim = TRUE,
                  longname = "Year")

scen <- ncdim_def(name = "SSP", 
                  units = "scenario", 
                  vals = c(1, 2, 3, 5), 
                  unlim = TRUE,
                  longname = "SSP scenario")

var_prec <- ncvar_def(name = "population",
                      units = "count",
                      dim = list(lon, lat, scen, time),
                      longname = "population",
                      missval = mv,
                      compression = 9)

ncout <- nc_create(filename, list(var_prec), force_v4 = TRUE)

ncatt_put(ncout, 0, "Title", "population")
ncatt_put(ncout, 0, "Source", "Gao et al. 2021")
ncatt_put(ncout, 0, "References", "Falchetta et al. (2023)")
ncatt_put(ncout, 0, "Created on", date())

for (i in c(1:5)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 1, match(i, c(1:5))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(6:10)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 2, match(i, c(6:10))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(11:15)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 3, match(i, c(11:15))), 
            count = c(-1, -1, 1, 1))
}

for (i in c(16:20)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = var_prec, 
            vals = values(shape_ely_diff[[i]]), 
            start = c(1, 1, 4, match(i, c(16:20))), 
            count = c(-1, -1, 1, 1))
}

nc_close(ncout)

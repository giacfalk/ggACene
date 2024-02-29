
line = .2
cex = 1
side = 3
adj=-0.065


par(mfrow = c(2, 2), mai = c(0.3, 0.25, 0.25, 0.15))


f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
f1_c = f1
f1_c[f1_c<0.01*100] <- NA
f1_c[f1 >= 0.01*100 & f1 < 0.1*100] <- 1
f1_c[f1 >= 0.1*100 & f1 < 0.25*100] <- 2
f1_c[f1 >= 0.25*100 & f1 < 0.5*100] <- 3
f1_c[f1 >= 0.5*100 & f1 < 0.75*100] <- 4
f1_c[f1 >= 0.75*100] <- 5
f1_c <- ratify(f1_c) # https://stackoverflow.com/questions/23840178/how-to-write-a-raster-with-rat-factors-in-r-raster-package
rat <- levels(f1_c)[[1]]#get the values of the unique cell frot the attribute table
rat$legend <- c("1-10%", "10-25%", "25-50%", "50-75%", ">75%")[1:nrow(levels(f1_c)[[1]])]
levels(f1_c) <- rat

# mask out pixels without population
pop_ssps_data <- stack(pop_ssps[2])[[10]]
f1_c <- raster::mask(f1_c, pop_ssps_data>1000, maskvalue=0)

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=c('white', "lightyellow", '#ffce61', '#93003a', 'black'), main="AC penetration (% of households), 2020" ,xaxt = "n", yaxt = "n", legend = FALSE, axes=FALSE, box=FALSE)
#legend(x='bottom', legend =rat$legend,fill = c('white', "lightyellow", '#ffce61', '#93003a', 'black'), horiz=T, box.col="transparent")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("A", side=side, line=line, cex=cex, adj=adj)

#

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
values(f1) <- ifelse(values(f1)<values(f1_bk), values(f1_bk), values(f1))

f1_c = f1
f1_c[f1_c<0.01*100] <- NA
f1_c[f1 >= 0.01*100 & f1 < 0.1*100] <- 1
f1_c[f1 >= 0.1*100 & f1 < 0.25*100] <- 2
f1_c[f1 >= 0.25*100 & f1 < 0.5*100] <- 3
f1_c[f1 >= 0.5*100 & f1 < 0.75*100] <- 4
f1_c[f1 >= 0.75*100] <- 5
f1_c <- ratify(f1_c) # https://stackoverflow.com/questions/23840178/how-to-write-a-raster-with-rat-factors-in-r-raster-package
rat <- levels(f1_c)[[1]]#get the values of the unique cell frot the attribute table
rat$legend <- c("1-10%", "10-25%", "25-50%", "50-75%", ">75%")[1:nrow(levels(f1_c)[[1]])]
levels(f1_c) <- rat

# mask out pixels without population
pop_ssps_data <- stack(pop_ssps[2])[[45]]
f1_c <- raster::mask(f1_c, pop_ssps_data>1000, maskvalue=0)
f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=c('white', "lightyellow", '#ffce61', '#93003a', 'black'), main="AC penetration (% of households), SSP245, 2050", xaxt = "n", yaxt = "n", legend = FALSE, axes=FALSE, box=FALSE)
legend(x='bottom', legend =rat$legend,fill = c('white', "lightyellow", '#ffce61', '#93003a', 'black'), horiz=T, box.col="transparent", text.width=42)
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("B", side=side, line=line, cex=cex, adj=adj)


f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
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

# mask out pixels without population
pop_ssps_data <- stack(pop_ssps[2])[[15]]
f1_c <- raster::mask(f1_c, pop_ssps_data>1000, maskvalue=0)

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=c('white', "lightyellow", '#ffce61', '#93003a', 'black'), main="AC electricity use (GWh/yr.), 2020", xaxt = "n", yaxt = "n", legend = FALSE, axes=FALSE, box=FALSE)
#legend(x='bottom', legend =rat$legend,fill = heat.colors(5, rev=T), horiz=T, box.col="transparent")
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("C", side=side, line=line, cex=cex, adj=adj)

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2050")
f1_bk = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
values(f1) <- ifelse(values(f1)<values(f1_bk), values(f1_bk), values(f1))
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

# mask out pixels without population
pop_ssps_data <- stack(pop_ssps[2])[[45]]
f1_c <- raster::mask(f1_c, pop_ssps_data>1000, maskvalue=0)

f1_c <- projectRaster(f1_c, crs=newproj)

plot(f1_c, col=c('white', "lightyellow", '#ffce61', '#93003a', 'black'), main="AC electricity use (GWh/yr.), SSP245, 2050", xaxt = "n", yaxt = "n", legend = FALSE, axes=FALSE, box=FALSE)
legend(x='bottom', legend =rat$legend,fill = c('white', "lightyellow", '#ffce61', '#93003a', 'black'), horiz=T, box.col="transparent", text.width = 42)
plot(wrld_simpl  %>% spTransform(newproj), add=TRUE, fill=NA, lwd=0.01)
mtext("D", side=side, line=line, cex=cex, adj=adj)

dev.off()

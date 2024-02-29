
library(rmapshaper)

pop_ssps <- list.files(path=paste0(wd, "supporting_data/pop_downscaled_spps"), recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- stack(pop_ssps[2])[[10]]

data(wrld_simpl)
wrld_simpl_sf <- st_as_sf(wrld_simpl)
wrld_simpl_sf <- st_transform(wrld_simpl_sf, "ESRI:54009")
wrld_simpl_sf <- filter(wrld_simpl_sf, NAME!="Antarctica")
wrld_simpl_sf <- ms_filter_islands(wrld_simpl_sf, min_area = 12391399903)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
f1[pop_ssps_data==0] <- NA

f1_c <- projectRaster(f1, crs=newproj)
  
f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

a <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
                theme_void() +
                theme(
                  legend.position = "bottom"
                )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name="%", limits=c(0, 100))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))

pop_ssps_data <- stack(pop_ssps[2])[[10]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100

f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X354.001344086022, na.rm=T))

b <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_x_continuous(limits=c(0, 75))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 100))+
  ylab("Latitude")+
  xlab("Avg. AC penetration (2020)")
  

###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

pop_ssps_data <- stack(pop_ssps[2])[[45]]
f1[pop_ssps_data==0] <- NA

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

c <-  ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name="%", limits=c(0, 100))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))

pop_ssps_data <- stack(pop_ssps[2])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2050")*100
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

d <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_x_continuous(limits=c(0, 75))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 100))+
  ylab("Latitude")+
  xlab("Avg. AC penetration (2050, SSP245)")

################

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

e <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name="GWh/yr.", trans="log10", limits=c(0.0001, 8514))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))


pop_ssps_data <- stack(pop_ssps[2])[[10]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X354.001344086022, na.rm=T))

f <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans="log10", limits=c(0.0001, 8514))+
  scale_x_continuous(limits=c(0, 3100))+
  ylab("Latitude")+
  xlab("Avg. AC electricity use (2020)")


###

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2050")
f1_bk = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

g <-  ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name="GWh/yr.", trans="log10", limits=c(0.0001, 8514))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))


pop_ssps_data <- stack(pop_ssps[2])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

h <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans="log10", limits=c(0.0001, 8514))+
  scale_x_continuous(limits=c(0, 3100))+
  ylab("Latitude")+
  xlab("Avg. AC electricity use (2050, SSP245)")


library(patchwork)

b + a + d + c + f + e + h + g + plot_layout(widths = rep(c(0.2, 1, 0.2, 1), 2), ncol = 4) + plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", "D", "")))

ggsave("results/graphs_tables/maps.pdf", height = 5*1.2, width = 5*2, scale=1.35)

#####
#####


library(rmapshaper)

pop_ssps <- list.files(path=paste0(wd, "supporting_data/pop_downscaled_spps"), recursive = T, pattern="nc", full.names = T)

data(wrld_simpl)
wrld_simpl_sf <- st_as_sf(wrld_simpl)
wrld_simpl_sf <- st_transform(wrld_simpl_sf, "ESRI:54009")
wrld_simpl_sf <- filter(wrld_simpl_sf, NAME!="Antarctica")
wrld_simpl_sf <- ms_filter_islands(wrld_simpl_sf, min_area = 12391399903)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP1.2050")*100
pop_ssps_data <- stack(pop_ssps[1])[[45]]
f1[pop_ssps_data==0] <- NA
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

a <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name="%", limits=c(0, 100))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))

pop_ssps_data <- stack(pop_ssps[1])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ac, pop_ssps_data, "SSP1.2050")*100
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

b <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_x_continuous(limits=c(0, 75))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 100))+
  ylab("Latitude")+
  xlab("Avg. AC penetration (2050, SSP126)")


###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP3.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP3.2020")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

pop_ssps_data <- stack(pop_ssps[3])[[45]]
f1[pop_ssps_data==0] <- NA
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

c <-  ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name="%", limits=c(0, 100))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))

pop_ssps_data <- stack(pop_ssps[3])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ac, pop_ssps_data, "SSP1.2050")*100
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

d <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_x_continuous(limits=c(0, 75))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 100))+
  ylab("Latitude")+
  xlab("Avg. AC penetration (2050, SSP370)")

###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP5.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP5.2020")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

pop_ssps_data <- stack(pop_ssps[5])[[45]]
f1[pop_ssps_data==0] <- NA
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

e <-  ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name="%", limits=c(0, 100))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))

pop_ssps_data <- stack(pop_ssps[5])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ac, pop_ssps_data, "SSP5.2050")*100
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

f <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_x_continuous(limits=c(0, 75))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits=c(0, 100))+
  ylab("Latitude")+
  xlab("Avg. AC penetration (2050, SSP585)")


################

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP1_2050")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

g <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name="GWh/yr.", trans="log10", limits=c(0.0001, 8514))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))


pop_ssps_data <- stack(pop_ssps[1])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP1_2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

h <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans="log10", limits=c(0.0001, 8514))+
  scale_x_continuous(limits=c(0, 3100))+
  ylab("Latitude")+
  xlab("Avg. AC electricity use (2050, SSP126)")


f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP3_2050")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

i <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name="GWh/yr.", trans="log10", limits=c(0.0001, 8514))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))


pop_ssps_data <- stack(pop_ssps[3])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP3_2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

j <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans="log10", limits=c(0.0001, 8514))+
  scale_x_continuous(limits=c(0, 3100))+
  ylab("Latitude")+
  xlab("Avg. AC electricity use (2050, SSP370)")


f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP5_2050")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)
f1_c <- na.omit(f1_c)

k <- ggplot()+
  geom_raster(data=f1_c, aes(x=x, y=y, fill=layer))+
  theme_void() +
  theme(
    legend.position = "bottom"
  )+
  geom_sf(data=wrld_simpl_sf, fill="transparent", colour="white", lwd=0.001)+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name="GWh/yr.", trans="log10", limits=c(0.0001, 8514))+ theme(
    legend.key.height = unit(1, "line"),  # Adjust the size of the legend keys
    legend.key.width = unit(2.5, "line"),   # Adjust the size of the legend keys
    legend.spacing.x = unit(0.15, "cm"),   # Adjust horizontal spacing between legend items
    legend.spacing.y = unit(0.15, "cm")      # Adjust vertical spacing between legend items
  )+
  coord_sf(xlim = c(-11960297, 16008503))


pop_ssps_data <- stack(pop_ssps[5])[[45]]
pop <- as.data.frame(pop_ssps_data, xy=T)
f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP5_2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(layer=weighted.mean(layer, X389.001344086022, na.rm=T))

l <- ggplot(aggdata) + 
  theme_classic()+
  geom_col(aes(x=layer, y=y, fill=layer), show.legend = F) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1, trans="log10", limits=c(0.0001, 8514))+
  scale_x_continuous(limits=c(0, 3100))+
  ylab("Latitude")+
  xlab("Avg. AC electricity use (2050, SSP585)")


library(patchwork)

b + a + h + g + d + c + j + i + f + e  +  l + k + plot_layout(widths = rep(c(0.2, 1, 0.2, 1), 3), ncol = 4) + plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "", "D", "", "E", "", "F", "")))

ggsave("results/graphs_tables/maps_si.pdf", height = 5*1.5, width = 5*2, scale=1.25)

#####
#####
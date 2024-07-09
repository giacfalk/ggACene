
# shape_ac <- readRDS("output_data/shape_ac.Rds")
# shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")
# shape_ely_ac<- readRDS("output_data/shape_ely_ac.Rds")
# shape_ely_noac<- readRDS("output_data/shape_ely_noac.Rds")

library(rmapshaper)

pop_ssps <- list.files(path=paste0(wd, "/supporting_data/pop_downscaled_spps"), recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- stack(pop_ssps[2])[[10]]

data(wrld_simpl)
wrld_simpl_sf <- st_as_sf(wrld_simpl)
wrld_simpl_sf <- st_transform(wrld_simpl_sf, "ESRI:54009")
wrld_simpl_sf <- filter(wrld_simpl_sf, NAME!="Antarctica")
wrld_simpl_sf <- ms_filter_islands(wrld_simpl_sf, min_area = 12391399903)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100

newproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X354.001344086022<1000] <- NA

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

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(X354.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X354.001344086022*(layer/100), na.rm=T)) %>% ungroup()

aggdata <- melt(aggdata, 1)

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

b <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e6, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  ylab("")+
  xlab("")+
  scale_x_continuous(limits = c(0, 175), position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(X354.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X354.001344086022*(layer/100), na.rm=T))

aggdata2 <- melt(aggdata2, 1)

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x<- as.factor(aggdata2$x)

b2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(y=value/1e6, x=x, fill=variable), show.legend = T) +
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  scale_x_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0, 75), position = "left")+
  theme(axis.text.y = element_text(size=11), aspect.ratio =.125)

###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[2])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

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

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T)) %>% ungroup()

aggdata <- melt(aggdata, 1)

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

d <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e6, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  ylab("Latitude")+
  xlab("")+
  scale_x_continuous(limits = c(0, 175), position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T))

aggdata2 <- melt(aggdata2, 1)

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x<- as.factor(aggdata2$x)

d2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(y=value/1e6, x=x, fill=variable), show.legend = T) +
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  scale_x_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0, 75), position = "left")+
  theme(axis.text.y = element_text(size=11), aspect.ratio =.125)

################

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2020")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[2])[[10]]

pop_ssps_data  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X354.001344086022<1000] <- NA


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


shape_ely_ac <- sf::st_set_geometry(shape_ely_ac, shape_ac$geometry)
shape_ely_noac <- sf::st_set_geometry(shape_ely_noac, shape_ac$geometry)

##

pop_ssps_data <- stack(pop_ssps[2])[[10]]
hhsize_r = fasterize(shape_ac, pop_ssps_data, "hhsize")
pop_ssps_data <- pop_ssps_data / hhsize_r

pop <- as.data.frame(pop_ssps_data, xy=T)
colnames(pop)[3] <- "X354.001344086022"

f1 = fasterize(shape_ely_ac, pop_ssps_data, "SSP2.2020")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- na.omit(f1_c)

f1 = fasterize(shape_ely_noac, pop_ssps_data, "SSP2.2020")
f1_c2 <- as.data.frame(f1, xy=T)
f1_c2 <- na.omit(f1_c2)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2020")*100
f1_c3 <- as.data.frame(f1, xy=T)
f1_c3 <- na.omit(f1_c3)

f1_c$layer2 <- f1_c2$layer
f1_c$layer3 <- f1_c3$layer

f1_c <- merge(f1_c, pop, by=c("x", "y"))

library(spatstat)

aggdata <- f1_c %>% filter(X354.001344086022>1000) %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(layer2 * X354.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X354.001344086022*(layer3/100), na.rm=T))

aggdata$X354.001344086022 <- NULL

aggdata <- melt(aggdata, 1)

aggdata$variable <- factor(aggdata$variable, levels=c("Without AC", "With AC"))

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)


f <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Latitude")+
  xlab("")+
  guides(fill = guide_legend(nrow = 1))+
  scale_x_continuous(position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)


aggdata2 <- f1_c %>% filter(X354.001344086022>1000) %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(layer2 * X354.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X354.001344086022*(layer3/100), na.rm=T))

aggdata2$X354.001344086022 <- NULL

aggdata2 <- melt(aggdata2, 1)

aggdata2$variable <- factor(aggdata2$variable, levels=c("Without AC", "With AC"))

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x <- as.factor(aggdata2$x)

f2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=x, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-180, 0, 180), labels=c(-180, 0, 180))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Longitude")+
  xlab("")+
  coord_flip()+
  scale_x_continuous(position = "bottom")+
  theme(axis.text.y = element_text(size=11), aspect.ratio = 0.125)

###

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2050")
f1_bk = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP2_2050")
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[2])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

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
hhsize_r = fasterize(shape_ac, pop_ssps_data, "hhsize")
pop_ssps_data <- pop_ssps_data / hhsize_r

pop <- as.data.frame(pop_ssps_data, xy=T)
colnames(pop)[3] <- "X389.001344086022"

f1 = fasterize(shape_ely_ac, pop_ssps_data, "SSP2.2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- na.omit(f1_c)

f1 = fasterize(shape_ely_noac, pop_ssps_data, "SSP2.2050")
f1_c2 <- as.data.frame(f1, xy=T)
f1_c2 <- na.omit(f1_c2)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP2.2050")*100
f1_c3 <- as.data.frame(f1, xy=T)
f1_c3 <- na.omit(f1_c3)

f1_c$layer2 <- f1_c2$layer
f1_c$layer3 <- f1_c3$layer

f1_c <- merge(f1_c, pop, by=c("x", "y"))

library(spatstat)

aggdata <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata$X389.001344086022 <- NULL

aggdata <- melt(aggdata, 1)

aggdata$variable <- factor(aggdata$variable, levels=c("Without AC", "With AC"))

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)


h <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Latitude")+
  xlab("")+
  guides(fill = guide_legend(nrow = 1))+
  scale_x_continuous(position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)


aggdata2 <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata2$X389.001344086022 <- NULL

aggdata2 <- melt(aggdata2, 1)

aggdata2$variable <- factor(aggdata2$variable, levels=c("Without AC", "With AC"))

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x <- as.factor(aggdata2$x)

h2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=x, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-180, 0, 180), labels=c(-180, 0, 180))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Longitude")+
  xlab("")+
  coord_flip()+
  scale_x_continuous(position = "bottom")+
  theme(axis.text.y = element_text(size=11), aspect.ratio = 0.125)

library(patchwork)

ggsave("results/graphs_tables/maps_1.pdf", a + theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_1a.pdf", b+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_1b.pdf", b2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_2.pdf", c+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_2a.pdf", d+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_2b.pdf", d2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_3.pdf", e+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_3a.pdf", f+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_3b.pdf", f2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_4.pdf", g+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_4a.pdf", h+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_4b.pdf", h2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/legends.pdf", cowplot::plot_grid(cowplot::get_legend(a), cowplot::get_legend(e), cowplot::get_legend(b), cowplot::get_legend(f), ncol = 2), height = 0.5, width = 4, scale=2.5)

#####
#####

# do the same for the SI!!!
pop_ssps <- list.files(path=paste0(wd, "/supporting_data/pop_downscaled_spps"), recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- stack(pop_ssps[1])[[45]]

data(wrld_simpl)
wrld_simpl_sf <- st_as_sf(wrld_simpl)
wrld_simpl_sf <- st_transform(wrld_simpl_sf, "ESRI:54009")
wrld_simpl_sf <- filter(wrld_simpl_sf, NAME!="Antarctica")
wrld_simpl_sf <- ms_filter_islands(wrld_simpl_sf, min_area = 12391399903)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP1.2050")*100

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

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

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T)) %>% ungroup()

aggdata <- melt(aggdata, 1)

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

b <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e6, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  ylab("")+
  xlab("")+
  scale_x_continuous(limits = c(0, 175), position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T))

aggdata2 <- melt(aggdata2, 1)

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x<- as.factor(aggdata2$x)

b2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(y=value/1e6, x=x, fill=variable), show.legend = T) +
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  scale_x_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0, 75), position = "left")+
  theme(axis.text.y = element_text(size=11), aspect.ratio =.125)


###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP3.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP3.2050")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[3])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

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
f1 = fasterize(shape_ac, pop_ssps_data, "SSP3.2050")*100
f1_c <- as.data.frame(f1, xy=T)
f1_c <- merge(f1_c, pop, by=c("x", "y"))
f1_c <- na.omit(f1_c)

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T)) %>% ungroup()

aggdata <- melt(aggdata, 1)

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

b <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e6, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  ylab("")+
  xlab("")+
  scale_x_continuous(limits = c(0, 175), position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T))

aggdata2 <- melt(aggdata2, 1)

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x<- as.factor(aggdata2$x)

d2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(y=value/1e6, x=x, fill=variable), show.legend = T) +
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  scale_x_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0, 75), position = "left")+
  theme(axis.text.y = element_text(size=11), aspect.ratio =.125)

###

f1 = fasterize(shape_ac, pop_ssps_data, "SSP5.2050")*100
f1_bk = fasterize(shape_ac, pop_ssps_data, "SSP5.2050")*100
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[5])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

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

aggdata <- f1_c %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T)) %>% ungroup()

aggdata <- melt(aggdata, 1)

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

f <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e6, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  ylab("")+
  xlab("")+
  scale_x_continuous(limits = c(0, 175), position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(X389.001344086022*((100-layer)/100), na.rm=T), `With AC`=sum(X389.001344086022*(layer/100), na.rm=T))

aggdata2 <- melt(aggdata2, 1)

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x<- as.factor(aggdata2$x)

f2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(y=value/1e6, x=x, fill=variable), show.legend = T) +
  scale_fill_manual(name="",values=c("darkred", "forestgreen"))+
  scale_x_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(0, 75), position = "left")+
  theme(axis.text.y = element_text(size=11), aspect.ratio =.125)

################

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP1_2050")
f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[1])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA


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
hhsize_r = fasterize(shape_ac, pop_ssps_data, "hhsize")
pop_ssps_data <- pop_ssps_data / hhsize_r

pop <- as.data.frame(pop_ssps_data, xy=T)
colnames(pop)[3] <- "X389.001344086022"

f1 = fasterize(shape_ely_ac, pop_ssps_data, "SSP1.2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- na.omit(f1_c)

f1 = fasterize(shape_ely_noac, pop_ssps_data, "SSP1.2050")
f1_c2 <- as.data.frame(f1, xy=T)
f1_c2 <- na.omit(f1_c2)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP1.2050")*100
f1_c3 <- as.data.frame(f1, xy=T)
f1_c3 <- na.omit(f1_c3)

f1_c$layer2 <- f1_c2$layer
f1_c$layer3 <- f1_c3$layer

f1_c <- merge(f1_c, pop, by=c("x", "y"))

library(spatstat)

aggdata <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata$X389.001344086022 <- NULL

aggdata <- melt(aggdata, 1)

aggdata$variable <- factor(aggdata$variable, levels=c("Without AC", "With AC"))

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)

h <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Latitude")+
  xlab("")+
  guides(fill = guide_legend(nrow = 1))+
  scale_x_continuous(position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)

aggdata2 <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata2$X354.001344086022 <- NULL

aggdata2 <- melt(aggdata2, 1)

aggdata2$variable <- factor(aggdata2$variable, levels=c("Without AC", "With AC"))

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x <- as.factor(aggdata2$x)

h2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=x, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-180, 0, 180), labels=c(-180, 0, 180))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Longitude")+
  xlab("")+
  coord_flip()+
  scale_x_continuous(position = "bottom")+
  theme(axis.text.y = element_text(size=11), aspect.ratio = 0.125)

###

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP3_2050")
f1_bk = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP3_2020")
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[3])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

i <-  ggplot()+
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
hhsize_r = fasterize(shape_ac, pop_ssps_data, "hhsize")
pop_ssps_data <- pop_ssps_data / hhsize_r

pop <- as.data.frame(pop_ssps_data, xy=T)
colnames(pop)[3] <- "X389.001344086022"

f1 = fasterize(shape_ely_ac, pop_ssps_data, "SSP3.2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- na.omit(f1_c)

f1 = fasterize(shape_ely_noac, pop_ssps_data, "SSP3.2050")
f1_c2 <- as.data.frame(f1, xy=T)
f1_c2 <- na.omit(f1_c2)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP3.2050")*100
f1_c3 <- as.data.frame(f1, xy=T)
f1_c3 <- na.omit(f1_c3)

f1_c$layer2 <- f1_c2$layer
f1_c$layer3 <- f1_c3$layer

f1_c <- merge(f1_c, pop, by=c("x", "y"))

library(spatstat)

aggdata <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata$X389.001344086022 <- NULL

aggdata <- melt(aggdata, 1)

aggdata$variable <- factor(aggdata$variable, levels=c("Without AC", "With AC"))

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)


j <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Latitude")+
  xlab("")+
  guides(fill = guide_legend(nrow = 1))+
  scale_x_continuous(position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)


aggdata2 <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata2$X354.001344086022 <- NULL

aggdata2 <- melt(aggdata2, 1)

aggdata2$variable <- factor(aggdata2$variable, levels=c("Without AC", "With AC"))

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x <- as.factor(aggdata2$x)

j2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=x, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-180, 0, 180), labels=c(-180, 0, 180))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Longitude")+
  xlab("")+
  coord_flip()+
  scale_x_continuous(position = "bottom")+
  theme(axis.text.y = element_text(size=11), aspect.ratio = 0.125)

###

f1 = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP5_2050")
f1_bk = fasterize(shape_ely_diff, pop_ssps_data, "ely_total_SSP5_2050")
# values(f1) <- ifelse(values(f1)<values(f1_bk) , values(f1_bk) , values(f1))

f1_c <- projectRaster(f1, crs=newproj)

f1_c <- as.data.frame(f1_c, xy=T)

pop_ssps_data <- stack(pop_ssps[5])[[45]]

pop_ssps_dataa  <- as.data.frame(projectRaster(pop_ssps_data, crs=newproj), xy=T)

f1_c <- na.omit(f1_c)

f1_c <- merge(f1_c, pop_ssps_dataa, by=c("x", "y"))

f1_c$layer[f1_c$X389.001344086022<1000] <- NA

k <-  ggplot()+
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
hhsize_r = fasterize(shape_ac, pop_ssps_data, "hhsize")
pop_ssps_data <- pop_ssps_data / hhsize_r

pop <- as.data.frame(pop_ssps_data, xy=T)
colnames(pop)[3] <- "X389.001344086022"

f1 = fasterize(shape_ely_ac, pop_ssps_data, "SSP5.2050")
f1_c <- as.data.frame(f1, xy=T)
f1_c <- na.omit(f1_c)

f1 = fasterize(shape_ely_noac, pop_ssps_data, "SSP5.2050")
f1_c2 <- as.data.frame(f1, xy=T)
f1_c2 <- na.omit(f1_c2)

f1 = fasterize(shape_ac, pop_ssps_data, "SSP5.2050")*100
f1_c3 <- as.data.frame(f1, xy=T)
f1_c3 <- na.omit(f1_c3)

f1_c$layer2 <- f1_c2$layer
f1_c$layer3 <- f1_c3$layer

f1_c <- merge(f1_c, pop, by=c("x", "y"))

library(spatstat)

aggdata <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(y) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata$X389.001344086022 <- NULL

aggdata <- melt(aggdata, 1)

aggdata$variable <- factor(aggdata$variable, levels=c("Without AC", "With AC"))

aggdata <- filter(aggdata, y>-50 & y<90)

aggdata$y <- as.factor(aggdata$y)


l <- ggplot(aggdata) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=y, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-50, 0, 50), labels=c(-50, 0, 50))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Latitude")+
  xlab("")+
  guides(fill = guide_legend(nrow = 1))+
  scale_x_continuous(position = "top")+
  theme(axis.text.x = element_text(size=8), aspect.ratio =3.5)


aggdata2 <- f1_c %>% filter(X389.001344086022>1000) %>% group_by(x) %>% dplyr::summarise(`Without AC`=sum(layer2 * X389.001344086022*((100-layer3)/100), na.rm=T), `With AC`=sum(layer * X389.001344086022*(layer3/100), na.rm=T))

aggdata2$X354.001344086022 <- NULL

aggdata2 <- melt(aggdata2, 1)

aggdata2$variable <- factor(aggdata2$variable, levels=c("Without AC", "With AC"))

aggdata2 <- filter(aggdata2, x>-125 & x<175)

aggdata2$x <- as.factor(aggdata2$x)

l2 <- ggplot(aggdata2) + 
  theme_void()+
  geom_col(aes(x=value/1e9, y=x, fill=variable), show.legend = T) +
  scale_y_discrete(breaks=c(-180, 0, 180), labels=c(-180, 0, 180))+
  scale_fill_manual(name="",values=c("lightblue", "navyblue"))+
  ylab("Longitude")+
  xlab("")+
  coord_flip()+
  scale_x_continuous(position = "bottom")+
  theme(axis.text.y = element_text(size=11), aspect.ratio = 0.125)

library(patchwork)

ggsave("results/graphs_tables/maps_1_si.pdf", a + theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_1a_si.pdf", b+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_1b_si.pdf", b2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_2_si.pdf", c+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_2a_si.pdf", d+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_2b_si.pdf", d2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_3_si.pdf", e+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_3a_si.pdf", f+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_3b_si.pdf", f2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_4_si.pdf", g+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_4a_si.pdf", h+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_4b_si.pdf", h2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_5_si.pdf", i+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_5a_si.pdf", j+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_5b_si.pdf", j2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/maps_6_si.pdf", k+ theme(legend.position = "none"), height = 2, width = 4, scale=2.5)
ggsave("results/graphs_tables/maps_6a_si.pdf", l+ theme(legend.position = "none"), height = 2, width = 1.5, scale=2.5)
ggsave("results/graphs_tables/maps_6b_si.pdf", l2+ theme(legend.position = "none"), height = 1.5, width = 4, scale=2.5)

ggsave("results/graphs_tables/legends_si.pdf", cowplot::plot_grid(cowplot::get_legend(a), cowplot::get_legend(e), cowplot::get_legend(b), cowplot::get_legend(f), ncol = 2), height = 0.5, width = 4, scale=2.5)

# validation plot
#######

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

setwd(wd)

#####

shape_ac <- readRDS("output_data/shape_ac.Rds")
#shape_ely_diff<- readRDS("shape_ely_diff.Rds")

shape_ac$geometry <- NULL
#shape_ely_diff$geometry <- NULL

#shape_ac <- merge(shape_ac, shape_pop, "id")

shape_ac_s <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_ssp2_2020, na.rm=T))

#

# Load country data
load("results/global_wgt_dmcf.RData")

global = reg_ely$data
global <- dplyr::select(global, country, ac, weight)
global <- na.omit(global)
global <- dplyr::group_by(global, country) %>% dplyr::summarise(ac=weighted.mean(as.numeric(as.character(ac)), weight, na.rm=T))
global$country <- countrycode::countrycode(global$country, 'country.name', 'iso3c')

shape_ac_s <- merge(shape_ac_s, global, by.x="ISO3", by.y="country")

# produce validation plot

library(ggpubr)
library(ggrepel)

a_plot <- ggplot(shape_ac_s)+
  geom_abline()+
  geom_point(aes(x=SSP2.2010, y=ac))+
  geom_text_repel(aes(x=SSP2.2010, y=ac, label=ISO3), size=3)+
  annotate("text", x=0.05, y=1, label= paste0("r = ", round(cor(shape_ac_s$SSP2.2010, shape_ac_s$ac), 2)))+
  xlab("Modeled AC penetration rate, aggregated grid-cell estimates")+
  ylab("AC penetration rate from aggregated survey training data")+
  scale_x_continuous(limits=c(0,1),labels=scales::label_percent())+
  scale_y_continuous(limits=c(0,1),labels=scales::label_percent())


###################################
###################################
###################################
###################################

gdata::keep(a_plot, wd, sure = T) # Removes all previously created variables
gc()                  # frees up memory resources

setwd(wd)

#####

shape_ac <- readRDS("output_data/shape_ac.Rds")
#shape_ely_diff<- readRDS("shape_ely_diff.Rds")

shape_ac$geometry <- NULL
#shape_ely_diff$geometry <- NULL

#shape_ac <- merge(shape_ac, shape_pop, "id")

shape_ac_s <- dplyr::group_by(shape_ac, ISO3) %>% dplyr::summarise(SSP2.2010=weighted.mean(SSP2.2010, pop_ssps_data_hist, na.rm=T))

#

setwd(wd)

validation_data <- readxl::read_xlsx("supporting_data/ac_statistics.xlsx")


merger <- merge(shape_ac_s, validation_data, by.x="ISO3", by.y="country")

sample_c <- c("JPN", "USA", "CHN", "MEX", "BRA", "IDN", "IND", "ARG", "DEU", "GHA", "NGA", "PAK")

merger$InSample = as.factor(ifelse(merger$ISO3 %in% sample_c, "Yes", "No"))

merger$ISO3 <- ifelse(merger$source=="Davis_Gertler", paste0(merger$ISO3, "_DG"), paste0(merger$ISO3, "_IEA"))

# produce validation plot

library(ggpubr)
library(ggrepel)
#https://stackoverflow.com/questions/60143052/how-to-add-r2-for-each-facet-of-ggplot-in-r

#merger <- filter(merger, source!="Davis_Gertler")

b_plot <- ggplot(merger)+
  geom_abline()+
  geom_point(aes(x=SSP2.2010, y=ac, colour=InSample))+
  geom_text_repel(aes(x=SSP2.2010, y=ac, label=ISO3), size=3)+
  annotate("text", x=0.05, y=1, label= paste0("r = ", round(cor(merger$SSP2.2010, merger$ac), 2)))+
  xlab("Modeled AC penetration rate, aggregated grid-cell estimates")+
  ylab("AC penetration rate from alternative national statistics")+
  scale_x_continuous(limits=c(0,1),labels=scales::label_percent())+
  scale_y_continuous(limits=c(0,1),labels=scales::label_percent())

library(patchwork)

a_plot + b_plot

ggsave("results/graphs_tables/validation_with_nat_stats.png", scale=1.5, width = 7, height = 4)  

setwd(wd)

merger$bias <- merger$SSP2.2010 - merger$ac
merger <- dplyr::select(merger, ISO3, bias)
merger$ISO3 <- substr(merger$ISO3, 1, 3)
write_rds(merger, "supporting_data/bias.Rds")

############


#######

rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(haven)
library(readstata13)
library(stringr)
library(tidyverse)
library(sandwich)
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
library(sf)
library(pROC)
library(pbapply)

setwd(wd)

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

shape_ac[,grep("pop_", colnames(shape_ac))] <- shape_ac[,grep("pop_", colnames(shape_ac))] / shape_ac$hhsize

###

shape_ac <- shape_ac %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))

paths <- dplyr::select(shape_ac, id, region, decile, starts_with("SSP"))
paths <- dplyr::select(paths, -contains("q"))

paths <- reshape2::melt(paths, c(1, 2, 3))
paths$ssp <- substr(paths$variable, 1, 4)
paths$year <- as.numeric(substr(paths$variable, 6, 9))

paths$ISO3 <- as.character(paths$region)
paths$decile <- as.factor(paths$decile)
paths <- filter(paths, year==2020 | year==2050)

paths$level <- ifelse(paths$ssp=="SSP2" & paths$year==2020, "2020", ifelse(paths$ssp=="SSP2" & paths$year==2050, "2050, SSP245",  ifelse(paths$ssp=="SSP5" & paths$year==2050, "2050, SSP585", ifelse(paths$ssp=="SSP1" & paths$year==2050, "2050, SSP126", ifelse(paths$ssp=="SSP3" & paths$year==2050, "2050, SSP370", NA)))))

paths <- filter(paths, !is.na(level))
paths$level <- as.factor(paths$level)

popp <- dplyr::select(shape_ac %>% ungroup(), region, decile, starts_with("pop"))
popp$POP2005 <- NULL
popp <- reshape2::melt(popp, c(1, 2))
popp$ssp <- substr(popp$variable, 15, 18)
popp <- filter(popp, ssp !="hist")
popp$ssp <- toupper(popp$ssp)
popp$year <- as.numeric(substr(popp$variable, 20, 23))

popp$ISO3 <- as.character(popp$region)
popp$decile <- as.factor(popp$decile)
popp <- filter(popp, year==2020 | year==2050)

popp$level <- ifelse(popp$ssp=="SSP2" & popp$year==2020, "2020", ifelse(popp$ssp=="SSP2" & popp$year==2050, "2050, SSP245",  ifelse(popp$ssp=="SSP5" & popp$year==2050, "2050, SSP585", ifelse(popp$ssp=="SSP1" & popp$year==2050, "2050, SSP126", ifelse(popp$ssp=="SSP3" & popp$year==2050, "2050, SSP370", NA)))))

popp <- filter(popp, !is.na(level))
popp$level <- as.factor(popp$level)

colnames(popp)[4] <- "pop"

paths$pop <- popp$pop

globbo <- paths

globbo$ISO3 <- " GLOBAL"

paths <- bind_rows(paths, globbo)

#

View(paths %>% group_by(ISO3, decile, level) %>% dplyr::summarise(value=weighted.mean(value, pop, na.rm=T)*100))

#

lines_a <- ggplot(paths %>% filter(level=="2020" | level=="2050, SSP245"))+
  theme_classic()+
  geom_boxplot(aes(x=decile, y=(value)*100, weight=pop, group=interaction(decile, level), fill=level), colour="lightgrey", lwd=0.01, outlier.colour = "transparent")+
  facet_wrap(vars(ISO3), ncol=4, scales="free_y")+
  scale_fill_manual(name="", values=c("lightblue", "navyblue"))+
  xlab("Income quintile")+
  ylab("AC penetration rate")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region. ")+
  scale_alpha_discrete(range=c(0.1, 1), guide = 'none') +
  theme(strip.background = element_blank()) 

#####

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

shape_ely_diff_m <- dplyr::select(shape_ely_diff, id, starts_with("SSP"))
shape_ely_diff_m$region <- NULL
colnames(shape_ely_diff_m) <- gsub("SSP", "cons_AC_SSP", colnames(shape_ely_diff_m))

shape_ely_diff_m <- merge(shape_ac,shape_ely_diff_m, "id")

shape_ely_diff_m$pop_ssps_data_ssp2_2020 <- shape_ely_diff_m$pop_ssps_data_ssp2_2020 * shape_ely_diff_m$SSP2.2020
shape_ely_diff_m$pop_ssps_data_ssp2_2050 <- shape_ely_diff_m$pop_ssps_data_ssp2_2050 * shape_ely_diff_m$SSP2.2050

shape_ely_diff <- shape_ely_diff_m %>% group_by(region) %>% mutate(decile=ntile(exp(GDP_SSP2_2010), 5))
shape_ely_diff$decile <- as.factor(shape_ely_diff$decile)

paths2 <- dplyr::select(shape_ely_diff, id, region, decile, starts_with("cons_AC_SSP"))
paths2 <- dplyr::select(paths2, -contains("q"))

paths2 <- reshape2::melt(paths2, c(1, 2, 3))
paths2$ssp <- substr(paths2$variable, 9, 12)
paths2$year <- as.numeric(substr(paths2$variable, 14, 17))

paths2$ISO3 <- as.character(paths2$region)
paths2$decile <- as.factor(paths2$decile)
paths2 <- filter(paths2, year==2020 | year==2050)

paths2$level <- ifelse(paths2$ssp=="SSP2" & paths2$year==2020, "2020", ifelse(paths2$ssp=="SSP2" & paths2$year==2050, "2050, SSP245",  ifelse(paths2$ssp=="SSP5" & paths2$year==2050, "2050, SSP585", ifelse(paths2$ssp=="SSP1" & paths2$year==2050, "2050, SSP126", ifelse(paths2$ssp=="SSP3" & paths2$year==2050, "2050, SSP370", NA)))))

paths2 <- filter(paths2, !is.na(level))
paths2$level <- as.factor(paths2$level)

popp <- dplyr::select(shape_ely_diff %>% ungroup(), region, decile, starts_with("pop"))
popp$POP2005 <- NULL
popp <- reshape2::melt(popp, c(1, 2))
popp$ssp <- substr(popp$variable, 15, 18)
popp <- filter(popp, ssp !="hist")
popp$ssp <- toupper(popp$ssp)
popp$year <- as.numeric(substr(popp$variable, 20, 23))

popp$ISO3 <- as.character(popp$region)
popp$decile <- as.factor(popp$decile)
popp <- filter(popp, year==2020 | year==2050)

popp$level <- ifelse(popp$ssp=="SSP2" & popp$year==2020, "2020", ifelse(popp$ssp=="SSP2" & popp$year==2050, "2050, SSP245",  ifelse(popp$ssp=="SSP5" & popp$year==2050, "2050, SSP585", ifelse(popp$ssp=="SSP1" & popp$year==2050, "2050, SSP126", ifelse(popp$ssp=="SSP3" & popp$year==2050, "2050, SSP370", NA)))))

popp <- filter(popp, !is.na(level))
popp$level <- as.factor(popp$level)

colnames(popp)[4] <- "pop"

paths2$pop <- popp$pop

globbo <- paths2

globbo$ISO3 <- " GLOBAL"

paths2 <- bind_rows(paths2, globbo)


#######

lines_b <- ggplot(paths2 %>% filter(level=="2020" | level=="2050, SSP245") %>% group_by(decile, level) %>% mutate(
  q1 = quantile(value, 0.25, na.rm=T),
  q3 = quantile(value, 0.75, na.rm=T),
  iqr = q3 - q1,
  upper_bound = q3 + 1.5 * iqr,
  lower_bound = q1 - 1.5 * iqr
) %>%
  filter(value >= lower_bound & value <= upper_bound) %>%
  # Optional: Remove the added columns (q1, q3, etc.)
  dplyr::select(-q1, -q3, -iqr, -upper_bound, -lower_bound))+
  theme_classic()+
  geom_boxplot(aes(x=decile, y=(value), weight=pop, group=interaction(decile, level), fill=level), colour="lightgrey", lwd=0.01, outlier.colour = "transparent")+
  facet_wrap(vars(ISO3), ncol=4, scales="free_y")+
  scale_fill_manual(name="", values=c("lightblue", "navyblue"))+
  xlab("Income quintile")+
  ylab("AC electricity consumption (kWh/HH./yr.)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="horizontal")+
  labs(caption = "Income quintiles are calculated within each region. ")+
  scale_alpha_discrete(range=c(0.1, 1), guide = 'none') +
  theme(strip.background = element_blank())
  

library(patchwork)

(lines_a + lines_b) + plot_layout(guides = "collect", ncol = 1) + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", legend.direction = "horizontal")  & guides(colour = guide_legend(nrow = 1))

ggsave("results/graphs_tables/quntiles_plot.pdf", scale=1.25, height = 8, width = 7)


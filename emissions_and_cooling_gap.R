
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

# quantify cooling gap 

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL

###

# calculate IQR of GCMs

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("126", "245", "370", "585")){
  for(y in seq(2015, 2050, 1)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & grepl("CDD", names(shape_ac))
    
    varname <- paste0("CDD_", x, "_", y)
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.5)))
    
  }}

future_acc$id <- NULL

shape_ac <- bind_cols(shape_ac, future_acc)

###

shape_ac_s <- shape_ac
shape_ac_s[is.na(shape_ac_s)] = 0
shape_ac_s$SSP2.2050_cg <- (1- shape_ac_s$SSP2.2050)*shape_ac_s$pop_ssps_data_ssp2_2050
shape_ac_s$SSP5.2050_cg <- (1- shape_ac_s$SSP5.2050)*shape_ac_s$pop_ssps_data_ssp5_2050
shape_ac_s$SSP1.2050_cg <- (1- shape_ac_s$SSP1.2050)*shape_ac_s$pop_ssps_data_ssp1_2050
shape_ac_s$SSP3.2050_cg <- (1- shape_ac_s$SSP3.2050)*shape_ac_s$pop_ssps_data_ssp3_2050
shape_ac_s$hist_cg <- (1- shape_ac_s$SSP2.2020)*shape_ac_s$pop_ssps_data_ssp2_2020

shape_ac_s <- dplyr::select(shape_ac_s, region, CDD_126_2020, CDD_126_2050, CDD_245_2020, CDD_245_2050, CDD_370_2020, CDD_370_2050, CDD_585_2050, SSP1.2050_cg, SSP2.2050_cg, SSP3.2050_cg, SSP5.2050_cg, hist_cg, ISO3)

shape_ac_s <- bind_rows(shape_ac_s, shape_ac_s %>% mutate(region=" GLOBAL"))

shape_ac_s_cs_1 <- shape_ac_s %>% group_by(region) %>%  dplyr::arrange(CDD_245_2020) %>% dplyr::mutate(hist_cg_cs = cumsum(coalesce(hist_cg, 0))/sum(hist_cg, na.rm=T))

shape_ac_s_cs_2 <- shape_ac_s %>% group_by(region) %>% dplyr::arrange(CDD_245_2050) %>% dplyr::mutate(SSP2.2050_cs = cumsum(coalesce(SSP2.2050_cg, 0))/sum(SSP2.2050_cg, na.rm=T))

shape_ac_s_cs_3 <- shape_ac_s %>% group_by(region) %>% dplyr::arrange(CDD_585_2050) %>% dplyr::mutate(SSP5.2050_cs = cumsum(coalesce(SSP5.2050_cg, 0))/sum(SSP5.2050_cg, na.rm=T))

shape_ac_s_cs_4 <- shape_ac_s %>% group_by(region) %>% dplyr::arrange(CDD_126_2050) %>% dplyr::mutate(SSP1.2050_cs = cumsum(coalesce(SSP1.2050_cg, 0))/sum(SSP1.2050_cg, na.rm=T))

shape_ac_s_cs_5 <- shape_ac_s %>% group_by(region) %>% dplyr::arrange(CDD_370_2050) %>% dplyr::mutate(SSP3.2050_cs = cumsum(coalesce(SSP3.2050_cg, 0))/sum(SSP3.2050_cg, na.rm=T))

shape_ac_s_cs_1 <- shape_ac_s_cs_1 %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T))

cg_a <- ggplot()+
  theme_classic()+
  geom_vline(data=shape_ac_s_cs_1, aes(xintercept = thresh), linetype="dashed", alpha=.5, colour="black", lwd=0.7)+
  geom_line(data=shape_ac_s_cs_1, aes(y=hist_cg_cs, x=CDD_245_2020), colour='grey', lwd=0.5)+
  geom_line(data=shape_ac_s_cs_4, aes(y=SSP1.2050_cs, x=CDD_126_2050), colour="#fcfc65", lwd=0.5)+
    geom_line(data=shape_ac_s_cs_2, aes(y=SSP2.2050_cs, x=CDD_245_2050), colour="#facf96", lwd=0.5)+
  geom_line(data=shape_ac_s_cs_5, aes(y=SSP3.2050_cs, x=CDD_370_2050), colour="#e38202", lwd=0.5)+
  geom_line(data=shape_ac_s_cs_3, aes(y=SSP5.2050_cs, x=CDD_585_2050), colour="#b51209", lwd=0.5)+
  scale_y_continuous(labels=scales::label_percent())+
  facet_wrap(vars(region), nrow=2) +
  theme(strip.background = element_blank())+
  labs(caption = "Vertical dashed lines represent region-specific average historical CDDs/yr.")+
  xlab("CDDs / yr.")+
  ylab("Cumulative % of pop. without AC")+
  scale_color_manual(name='Scenario',
                     breaks=c('2020', 'SSP126', 'SSP245','SSP370', 'SSP585'),
                     values=c('2020'='grey', 'SSP126' = '#fcfc65', 'SSP245'='#facf96', 'SSP370' = '#e38202', 'SSP585'='#7d0404'))#

shape_ac_s_cs_1 <- shape_ac_s %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T))%>% dplyr::summarise(cg = sum(hist_cg[CDD_245_2020>thresh], na.rm=T))

shape_ac_s_cs_2 <- shape_ac_s %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T)) %>% dplyr::summarise(threshold=mean(thresh, na.rm=T), cg = sum(SSP2.2050_cg[CDD_245_2050>thresh], na.rm=T))

shape_ac_s_cs_3 <- shape_ac_s %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T))%>% dplyr::summarise(cg = sum(SSP5.2050_cg[CDD_585_2050>thresh], na.rm=T))

shape_ac_s_cs_4 <- shape_ac_s %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T)) %>% dplyr::summarise(threshold=mean(thresh, na.rm=T), cg = sum(SSP1.2050_cg[CDD_126_2050>thresh], na.rm=T))

shape_ac_s_cs_5 <- shape_ac_s %>% dplyr::group_by(region) %>% dplyr::mutate(thresh=weighted.mean(CDD_245_2020, hist_cg, na.rm=T))%>% dplyr::summarise(cg = sum(SSP3.2050_cg[CDD_370_2050>thresh], na.rm=T))

shape_ac_s_cs_1$scenario <- "2020"
shape_ac_s_cs_2$scenario <- "SSP245"
shape_ac_s_cs_3$scenario <- "SSP585"
shape_ac_s_cs_4$scenario <- "SSP126"
shape_ac_s_cs_5$scenario <- "SSP370"

shape_ac_s_cs <- bind_rows(shape_ac_s_cs_1, shape_ac_s_cs_2, shape_ac_s_cs_3, shape_ac_s_cs_4, shape_ac_s_cs_5)

shape_ac_s_cs <- as.data.frame(shape_ac_s_cs)

shape_ac_s_cs$scenario <- factor(shape_ac_s_cs$scenario, levels = c("2020", "SSP126", "SSP245", "SSP370", "SSP585"))

cg_b<- ggplot()+
  theme_classic()+
  geom_col(data=shape_ac_s_cs, aes(y=cg/1e9, x=region, fill=scenario), colour="black", position = "dodge")+
  scale_y_continuous("Billion people without AC \nexposed to CDDs/yr.> regional avg.)")+
  scale_fill_manual(name='Scenario',
                     breaks=c("2020", "SSP126", "SSP245", "SSP370", "SSP585"),
                     values=c('2020'='grey', 'SSP126' = '#fcfc65', 'SSP245'='#facf96', 'SSP370' = '#e38202', 'SSP585'='#7d0404'))+
  xlab("")+
  theme(legend.position = "bottom", legend.direction = "horizontal", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

library(patchwork)

cg_a + cg_b + plot_layout(ncol=1, heights = c(1, 1)) + plot_annotation(tag_levels = "A")

ggsave("results/graphs_tables/cooling_gap.pdf", scale=1.5, height = 6, width = 5)

#################

shape_ely_diff <- shape_ely_diff %>% dplyr::select(id, contains(".20") & contains("SSP"))
shape_ely_diff$region <- NULL
colnames(shape_ely_diff)[2:261] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:261])

shape_ely_diff <- ungroup(shape_ely_diff)

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ely_diff <- bind_cols(shape_ely_diff, future_acc_q1, future_acc_q3)

######

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ac$geometry <- NULL

## hhsize

shape_ac[,grep("pop_", colnames(shape_ac))] <- shape_ac[,grep("pop_", colnames(shape_ac))] / shape_ac$hhsize

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ac <- bind_cols(shape_ac, future_acc_q1, future_acc_q3)

#################

# emissions figure

require(data.table)
ci <- fread(paste0(wd, "/supporting_data/AR6_Scenarios_Database_ISO3_v1.0.csv"), header = T)

emis_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][88]
ci_emis <- filter(ci, Variable %in% emis_lab)
ci_emis <- filter(ci_emis, Scenario %in% c("SSP2_BASE", "SSP5-baseline", "SSP1-baseline", "SSP3-baseline"))
ci_emis <- dplyr::select(ci_emis, Scenario, Region, Unit, `2020`, `2050`)
ci_emis$`2050` <- ci_emis$`2050` * 1000000000000  # convert to grams
ci_emis$`2020` <- ci_emis$`2020` * 1000000000000  # convert to grams

ely_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][136:152]
ci_ely <- filter(ci, Variable %in% ely_lab[c(1:5, 9:10, 12, 15:17)])
ci_ely <- filter(ci_ely, Scenario %in% c("SSP2_BASE", "SSP5-baseline", "SSP1-baseline", "SSP3-baseline"))
ci_ely <- group_by(ci_ely, Scenario, Region, Unit) %>% dplyr::summarise(`2020` = sum(`2020`, na.rm=T), `2050`=sum(`2050`, na.rm=T))
ci_ely$`2050` <- ci_ely$`2050` * 277777777777.78 # convert to kWh
ci_ely$`2020` <- ci_ely$`2020` * 277777777777.78 # convert to kWh

ci_emis <- ci_emis[with(ci_emis, order(Region, Scenario)), ]
ci_ely <- ci_ely[with(ci_ely, order(Region, Scenario)), ]

ci_emis$ci_intensity_2020 <- ci_emis$`2020` / ci_ely$`2020`
ci_emis$ci_intensity_2050 <- ci_emis$`2050` / ci_ely$`2050`

# parse regions

shape_ely_diff2 <- shape_ely_diff %>% dplyr::select(1, 242:301)
shape_ely_diff_m <- merge(shape_ac,shape_ely_diff2, "id")

co2_glob1 <- dplyr::group_by(shape_ely_diff_m, region) %>% dplyr::summarise(cons_AC_SSP2.2010=sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=sum(cons_AC_SSP2.2030*SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=sum(cons_AC_SSP2.2040*SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=sum(cons_AC_SSP5.2010*SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=sum(cons_AC_SSP5.2020*SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=sum(cons_AC_SSP5.2030*SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=sum(cons_AC_SSP5.2040*SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q1=sum(cons_AC_SSP2.2010_q1*SSP2.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q1=sum(cons_AC_SSP2.2020_q1*SSP2.2020_q1*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q1=sum(cons_AC_SSP2.2030_q1*SSP2.2030_q1*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q1=sum(cons_AC_SSP2.2040_q1*SSP2.2040_q1*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q1=sum(cons_AC_SSP2.2050_q1*SSP2.2050_q1*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q1=sum(cons_AC_SSP5.2010_q1*SSP5.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q1=sum(cons_AC_SSP5.2020_q1*SSP5.2020_q1*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q1=sum(cons_AC_SSP5.2030_q1*SSP5.2030_q1*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q1=sum(cons_AC_SSP5.2040_q1*SSP5.2040_q1*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q1=sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q3=sum(cons_AC_SSP2.2010_q3*SSP2.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q3=sum(cons_AC_SSP2.2020_q3*SSP2.2020_q3*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q3=sum(cons_AC_SSP2.2030_q3*SSP2.2030_q3*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q3=sum(cons_AC_SSP2.2040_q3*SSP2.2040_q3*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q3=sum(cons_AC_SSP2.2050_q3*SSP2.2050_q3*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q3=sum(cons_AC_SSP5.2010_q3*SSP5.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q3=sum(cons_AC_SSP5.2020_q3*SSP5.2020_q3*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q3=sum(cons_AC_SSP5.2030_q3*SSP5.2030_q3*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q3=sum(cons_AC_SSP5.2040_q3*SSP5.2040_q3*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q3=sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm=T))


co2_glob2 <- dplyr::group_by(shape_ely_diff_m, region) %>% dplyr::summarise(cons_AC_SSP1.2010=sum(cons_AC_SSP1.2010*SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=sum(cons_AC_SSP1.2020*SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=sum(cons_AC_SSP1.2030*SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=sum(cons_AC_SSP1.2040*SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=sum(cons_AC_SSP3.2010*SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=sum(cons_AC_SSP3.2020*SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=sum(cons_AC_SSP3.2030*SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=sum(cons_AC_SSP3.2040*SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q1=sum(cons_AC_SSP1.2010_q1*SSP1.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q1=sum(cons_AC_SSP1.2020_q1*SSP1.2020_q1*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q1=sum(cons_AC_SSP1.2030_q1*SSP1.2030_q1*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q1=sum(cons_AC_SSP1.2040_q1*SSP1.2040_q1*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q1=sum(cons_AC_SSP1.2050_q1*SSP1.2050_q1*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q1=sum(cons_AC_SSP3.2010_q1*SSP3.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q1=sum(cons_AC_SSP3.2020_q1*SSP3.2020_q1*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q1=sum(cons_AC_SSP3.2030_q1*SSP3.2030_q1*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q1=sum(cons_AC_SSP3.2040_q1*SSP3.2040_q1*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q1=sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q3=sum(cons_AC_SSP1.2010_q3*SSP1.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q3=sum(cons_AC_SSP1.2020_q3*SSP1.2020_q3*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q3=sum(cons_AC_SSP1.2030_q3*SSP1.2030_q3*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q3=sum(cons_AC_SSP1.2040_q3*SSP1.2040_q3*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q3=sum(cons_AC_SSP1.2050_q3*SSP1.2050_q3*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q3=sum(cons_AC_SSP3.2010_q3*SSP3.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q3=sum(cons_AC_SSP3.2020_q3*SSP3.2020_q3*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q3=sum(cons_AC_SSP3.2030_q3*SSP3.2030_q3*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q3=sum(cons_AC_SSP3.2040_q3*SSP3.2040_q3*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q3=sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm=T))

co2_glob <- merge(co2_glob1, co2_glob2, "region")

# merge

ci_emis <- pivot_longer(ci_emis, 6:7)
ci_emis$name <- gsub("ci_intensity_", "", ci_emis$name)
ci_emis$Scenario <- ifelse(ci_emis$Scenario=="SSP5-baseline", "SSP585 (2050)", ifelse(ci_emis$Scenario=="SSP2_BASE", "SSP245 (2050)", ifelse(ci_emis$Scenario=="SSP1-baseline", "SSP126 (2050)", "SSP370 (2050)")))
ci_emis$Scenario <- ifelse(ci_emis$name=="2020", "2020", ci_emis$Scenario)

colnames(ci_emis)[7] <- "emission_factor"

ci_emis <- dplyr::select(ci_emis, Scenario, Region, emission_factor)

co2_glob <- melt(co2_glob, 1)

co2_glob$Scenario <- ifelse(grepl("SSP2.2020", co2_glob$variable), "2020", NA) 
co2_glob$Scenario <- ifelse(grepl("SSP2.2050", co2_glob$variable), "SSP245 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP5.2050", co2_glob$variable), "SSP585 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP1.2050", co2_glob$variable), "SSP126 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP3.2050", co2_glob$variable), "SSP370 (2050)", co2_glob$Scenario ) 

co2_glob <- filter(co2_glob, !is.na(Scenario))

co2_glob$q <- NA
co2_glob$q <- ifelse(grepl("q1", co2_glob$variable), "Q1", co2_glob$q ) 
co2_glob$q <- ifelse(grepl("q3", co2_glob$variable), "Q3", co2_glob$q ) 
co2_glob$q <- ifelse(is.na(co2_glob$q), "Median", co2_glob$q ) 

co2_glob$variable <- NULL

co2_glob$Region <- ifelse(co2_glob$region=="East Asia & Pacific", "CHN", NA)
co2_glob$Region <- ifelse(co2_glob$region=="Europe & Central Asia", "EU", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Latin America & Caribbean", "BRA", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Middle East & North Africa", "TUR", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="North America", "USA", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="South Asia", "IND", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Sub-Saharan Africa", "ZAF", co2_glob$Region)

#

co2_glob <- pivot_wider(co2_glob, names_from = q, values_from = value)

#

merger <- merge(co2_glob, ci_emis, by=c("Scenario", "Region"), all.x=T)

merger <- group_by(merger,region) %>% dplyr::mutate(emission_factor = ifelse(is.na(emission_factor), mean(emission_factor, na.rm=T), emission_factor)) 

merger$mtco2_q5 <- merger$Median * merger$emission_factor / 1000000000000 # convert to Mt
merger$mtco2_q1 <- merger$Q1 * merger$emission_factor / 1000000000000 # convert to Mt
merger$mtco2_q3 <- merger$Q3 * merger$emission_factor / 1000000000000 # convert to Mt

merger <- dplyr::select(merger, Scenario, region, mtco2_q5, mtco2_q1, mtco2_q3)

#

merger <- merger %>% group_by(Scenario, region) %>% dplyr::summarise(mtco2_q5 = mean(mtco2_q5), mtco2_q1 = mean(mtco2_q1), mtco2_q3 = mean(mtco2_q3))

mergers <- merger

merger <- merger %>%
  bind_rows(data.frame(region = "Total", Scenario=unique(merger$Scenario), mtco2_q5 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q5=sum(mtco2_q5)) %>% pull(mtco2_q5), mtco2_q1 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q1=sum(mtco2_q1)) %>% pull(mtco2_q1), mtco2_q3 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q3=sum(mtco2_q3)) %>% pull(mtco2_q3)))
           
#

colnames(merger) <- c("Scenario", "Region", "Mt CO2", "Mt CO2 (Q1)", "Mt CO2 (Q3)")

merger$`Mt CO2` <- round(merger$`Mt CO2`, 1)
merger$`Mt CO2 (Q1)` <- round(merger$`Mt CO2 (Q1)`, 1)
merger$`Mt CO2 (Q3)` <- round(merger$`Mt CO2 (Q3)`, 1)

merger_bk <- merger

###

merger$`Mt CO2` <- paste0(merger$`Mt CO2`, " (", merger$`Mt CO2 (Q1)`, " - ", merger$`Mt CO2 (Q3)`, ")")

merger <- dplyr::select(merger, 1,2,3)

merger <- group_by(merger, Scenario, Region) %>% sample_n(1)

merger <- pivot_wider(merger, names_from = Scenario, values_from = `Mt CO2`)

sink("results/graphs_tables/emissions.tex")
xtable(merger)
sink()

###

merger_p <- merger_bk
merger_p$year <- ifelse(merger_p$Scenario=="2020", 2020, 2050)
merger_p$Scenario <- substr(merger_p$Scenario, 1, 6)
merger_p <- filter(merger_p, Region!="Total")

merger_p_glob <- merger_p %>% group_by(year, Scenario) %>% dplyr::summarise_if(is.numeric, sum, na.rm=T)

merger_p$type ="Region"
merger_p_glob$type ="Global"

merger_p_glob$Region=" GLOBAL"
merger_p <- bind_rows(merger_p, merger_p_glob)

colnames(merger_p)[4:5] <- c("mtco2_q1", "mtco2_q3")

emissions_figure <- ggplot(merger_p)+
  theme_classic()+
  ggtitle("Projected evolution of GHG emissions from residential AC")+
  geom_errorbar(aes(x=as.factor(year), ymin=mtco2_q1, ymax= mtco2_q3, group=Scenario, colour=Scenario), position = "dodge", width =0.2)+
  facet_wrap(vars(Region), scales = "free", nrow=2)+
  xlab("Year")+
  ylab("AC GHG emissions (Mt CO2e/yr.)")+
  labs(caption = "Ribbon: IQR of CMIP6 GCMs") +
  theme(strip.background = element_blank())+
  scale_colour_manual(name="Scenario", values=c("grey", "#fcfc65", "#facf96", "#e38202", "#7d0404"))

saveRDS(emissions_figure, "results/emissions_figure.rds")

################################################

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ely_diff<- readRDS("output_data/shape_ely_diff.Rds")

shape_ac$geometry <- NULL
shape_ely_diff$geometry <- NULL


shape_ely_diff <- shape_ely_diff %>% dplyr::select(id, contains(".20") & contains("SSP"))
shape_ely_diff$region <- NULL
colnames(shape_ely_diff)[2:261] <- paste0("cons_AC_", colnames(shape_ely_diff)[2:261])

shape_ely_diff <- ungroup(shape_ely_diff)

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ely_diff))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ely_diff)) & grepl(x, names(shape_ely_diff)) & (names(shape_ely_diff)!=paste0("cons_AC_", x, ".", y))
    
    varname <- paste0("cons_AC_", x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ely_diff %>% dplyr::select(colnames(shape_ely_diff)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ely_diff <- bind_cols(shape_ely_diff, future_acc_q1, future_acc_q3)

######

shape_ac <- readRDS("output_data/shape_ac.Rds")
shape_ac$geometry <- NULL

## hhsize

shape_ac[,grep("pop_", colnames(shape_ac))] <- shape_ac[,grep("pop_", colnames(shape_ac))] / shape_ac$hhsize

rowwise_quantile <- function(data, probs, na.rm=T) {
  apply(data, 1, quantile, probs = probs, na.rm=T)
}

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q1")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.25)))
    
  }}


future_acc$id <- NULL

future_acc_q1 <- future_acc

###

future_acc <- data.frame(id=1:nrow(shape_ac))

for(x in c("SSP1", "SSP2", "SSP3", "SSP5")){
  for(y in seq(2010, 2050, 10)){
    
    col1_columns <- grepl(y, names(shape_ac)) & grepl(x, names(shape_ac)) & (names(shape_ac)!=paste0(x, ".", y)) & (!grepl("GDP", names(shape_ac)) & !grepl("URB", names(shape_ac)) & !grepl("women", names(shape_ac)) & !grepl("age", names(shape_ac)) & !grepl("edu", names(shape_ac)))
    
    varname <- paste0(x, ".", y, "_q3")
    
    future_acc <-bind_cols(future_acc, shape_ac %>% dplyr::select(colnames(shape_ac)[col1_columns]) %>% dplyr::summarise(!!varname := rowwise_quantile(., 0.75)))
    
  }}


future_acc$id <- NULL
future_acc_q3 <- future_acc

shape_ac <- bind_cols(shape_ac, future_acc_q1, future_acc_q3)

require(data.table)
ci <- fread(paste0(wd, "/supporting_data/AR6_Scenarios_Database_ISO3_v1.0.csv"), header = T)

emis_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][88]
ci_emis <- filter(ci, Variable %in% emis_lab)
ci_emis <- filter(ci_emis, Scenario %in% c("SSP2_BASE", "SSP5-baseline", "SSP1-baseline", "SSP3-baseline"))
ci_emis <- dplyr::select(ci_emis, Scenario, Region, Unit, `2020`, `2050`)
ci_emis$`2050` <- ci_emis$`2050` * 1000000000000  # convert to grams
ci_emis$`2020` <- ci_emis$`2020` * 1000000000000  # convert to grams

ely_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][136:152]
ci_ely <- filter(ci, Variable %in% ely_lab[c(1:5, 9:10, 12, 15:17)])
ci_ely <- filter(ci_ely, Scenario %in% c("SSP2_BASE", "SSP5-baseline", "SSP1-baseline", "SSP3-baseline"))
ci_ely <- group_by(ci_ely, Scenario, Region, Unit) %>% dplyr::summarise(`2020` = sum(`2020`, na.rm=T), `2050`=sum(`2050`, na.rm=T))
ci_ely$`2050` <- ci_ely$`2050` * 277777777777.78 # convert to kWh
ci_ely$`2020` <- ci_ely$`2020` * 277777777777.78 # convert to kWh

ci_emis <- ci_emis[with(ci_emis, order(Region, Scenario)), ]
ci_ely <- ci_ely[with(ci_ely, order(Region, Scenario)), ]

ci_emis$ci_intensity_2020 <- ci_emis$`2020` / ci_ely$`2020`
ci_emis$ci_intensity_2050 <- ci_emis$`2050` / ci_ely$`2050`

# parse regions

shape_ely_diff2 <- shape_ely_diff %>% dplyr::select(1, 242:301)
shape_ely_diff_m <- merge(shape_ac,shape_ely_diff2, "id")

co2_glob1 <- dplyr::group_by(shape_ely_diff_m, region, ISO3) %>% dplyr::summarise(cons_AC_SSP2.2010=sum(cons_AC_SSP2.2010*SSP2.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020=sum(cons_AC_SSP2.2020*SSP2.2020*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030=sum(cons_AC_SSP2.2030*SSP2.2030*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040=sum(cons_AC_SSP2.2040*SSP2.2040*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050=sum(cons_AC_SSP2.2050*SSP2.2050*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010=sum(cons_AC_SSP5.2010*SSP5.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020=sum(cons_AC_SSP5.2020*SSP5.2020*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030=sum(cons_AC_SSP5.2030*SSP5.2030*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040=sum(cons_AC_SSP5.2040*SSP5.2040*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050=sum(cons_AC_SSP5.2050*SSP5.2050*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q1=sum(cons_AC_SSP2.2010_q1*SSP2.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q1=sum(cons_AC_SSP2.2020_q1*SSP2.2020_q1*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q1=sum(cons_AC_SSP2.2030_q1*SSP2.2030_q1*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q1=sum(cons_AC_SSP2.2040_q1*SSP2.2040_q1*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q1=sum(cons_AC_SSP2.2050_q1*SSP2.2050_q1*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q1=sum(cons_AC_SSP5.2010_q1*SSP5.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q1=sum(cons_AC_SSP5.2020_q1*SSP5.2020_q1*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q1=sum(cons_AC_SSP5.2030_q1*SSP5.2030_q1*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q1=sum(cons_AC_SSP5.2040_q1*SSP5.2040_q1*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q1=sum(cons_AC_SSP5.2050_q1*SSP5.2050_q1*(pop_ssps_data_ssp5_2050), na.rm=T), cons_AC_SSP2.2010_q3=sum(cons_AC_SSP2.2010_q3*SSP2.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP2.2020_q3=sum(cons_AC_SSP2.2020_q3*SSP2.2020_q3*(pop_ssps_data_ssp2_2020), na.rm=T), cons_AC_SSP2.2030_q3=sum(cons_AC_SSP2.2030_q3*SSP2.2030_q3*(pop_ssps_data_ssp2_2030), na.rm=T), cons_AC_SSP2.2040_q3=sum(cons_AC_SSP2.2040_q3*SSP2.2040_q3*(pop_ssps_data_ssp2_2040), na.rm=T), cons_AC_SSP2.2050_q3=sum(cons_AC_SSP2.2050_q3*SSP2.2050_q3*(pop_ssps_data_ssp2_2050), na.rm=T), cons_AC_SSP5.2010_q3=sum(cons_AC_SSP5.2010_q3*SSP5.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP5.2020_q3=sum(cons_AC_SSP5.2020_q3*SSP5.2020_q3*(pop_ssps_data_ssp5_2020), na.rm=T), cons_AC_SSP5.2030_q3=sum(cons_AC_SSP5.2030_q3*SSP5.2030_q3*(pop_ssps_data_ssp5_2030), na.rm=T), cons_AC_SSP5.2040_q3=sum(cons_AC_SSP5.2040_q3*SSP5.2040_q3*(pop_ssps_data_ssp5_2040), na.rm=T), cons_AC_SSP5.2050_q3=sum(cons_AC_SSP5.2050_q3*SSP5.2050_q3*(pop_ssps_data_ssp5_2050), na.rm=T))


co2_glob2 <- dplyr::group_by(shape_ely_diff_m, region, ISO3) %>% dplyr::summarise(cons_AC_SSP1.2010=sum(cons_AC_SSP1.2010*SSP1.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020=sum(cons_AC_SSP1.2020*SSP1.2020*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030=sum(cons_AC_SSP1.2030*SSP1.2030*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040=sum(cons_AC_SSP1.2040*SSP1.2040*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050=sum(cons_AC_SSP1.2050*SSP1.2050*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010=sum(cons_AC_SSP3.2010*SSP3.2010*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020=sum(cons_AC_SSP3.2020*SSP3.2020*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030=sum(cons_AC_SSP3.2030*SSP3.2030*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040=sum(cons_AC_SSP3.2040*SSP3.2040*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050=sum(cons_AC_SSP3.2050*SSP3.2050*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q1=sum(cons_AC_SSP1.2010_q1*SSP1.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q1=sum(cons_AC_SSP1.2020_q1*SSP1.2020_q1*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q1=sum(cons_AC_SSP1.2030_q1*SSP1.2030_q1*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q1=sum(cons_AC_SSP1.2040_q1*SSP1.2040_q1*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q1=sum(cons_AC_SSP1.2050_q1*SSP1.2050_q1*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q1=sum(cons_AC_SSP3.2010_q1*SSP3.2010_q1*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q1=sum(cons_AC_SSP3.2020_q1*SSP3.2020_q1*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q1=sum(cons_AC_SSP3.2030_q1*SSP3.2030_q1*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q1=sum(cons_AC_SSP3.2040_q1*SSP3.2040_q1*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q1=sum(cons_AC_SSP3.2050_q1*SSP3.2050_q1*(pop_ssps_data_ssp3_2050), na.rm=T), cons_AC_SSP1.2010_q3=sum(cons_AC_SSP1.2010_q3*SSP1.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP1.2020_q3=sum(cons_AC_SSP1.2020_q3*SSP1.2020_q3*(pop_ssps_data_ssp1_2020), na.rm=T), cons_AC_SSP1.2030_q3=sum(cons_AC_SSP1.2030_q3*SSP1.2030_q3*(pop_ssps_data_ssp1_2030), na.rm=T), cons_AC_SSP1.2040_q3=sum(cons_AC_SSP1.2040_q3*SSP1.2040_q3*(pop_ssps_data_ssp1_2040), na.rm=T), cons_AC_SSP1.2050_q3=sum(cons_AC_SSP1.2050_q3*SSP1.2050_q3*(pop_ssps_data_ssp1_2050), na.rm=T), cons_AC_SSP3.2010_q3=sum(cons_AC_SSP3.2010_q3*SSP3.2010_q3*(pop_ssps_data_hist), na.rm=T), cons_AC_SSP3.2020_q3=sum(cons_AC_SSP3.2020_q3*SSP3.2020_q3*(pop_ssps_data_ssp3_2020), na.rm=T), cons_AC_SSP3.2030_q3=sum(cons_AC_SSP3.2030_q3*SSP3.2030_q3*(pop_ssps_data_ssp3_2030), na.rm=T), cons_AC_SSP3.2040_q3=sum(cons_AC_SSP3.2040_q3*SSP3.2040_q3*(pop_ssps_data_ssp3_2040), na.rm=T), cons_AC_SSP3.2050_q3=sum(cons_AC_SSP3.2050_q3*SSP3.2050_q3*(pop_ssps_data_ssp3_2050), na.rm=T))

co2_glob <- merge(co2_glob1, co2_glob2, by=c("ISO3", "region"))

# merge

ci_emis <- pivot_longer(ci_emis, 6:7)
ci_emis$name <- gsub("ci_intensity_", "", ci_emis$name)
ci_emis$Scenario <- ifelse(ci_emis$Scenario=="SSP5-baseline", "SSP585 (2050)", ifelse(ci_emis$Scenario=="SSP2_BASE", "SSP245 (2050)", ifelse(ci_emis$Scenario=="SSP1-baseline", "SSP126 (2050)", "SSP370 (2050)")))
ci_emis$Scenario <- ifelse(ci_emis$name=="2020", "2020", ci_emis$Scenario)

colnames(ci_emis)[7] <- "emission_factor"

ci_emis <- dplyr::select(ci_emis, Scenario, Region, emission_factor)

co2_glob <- melt(co2_glob, 1:2)

co2_glob$Scenario <- ifelse(grepl("SSP2.2020", co2_glob$variable), "2020", NA) 
co2_glob$Scenario <- ifelse(grepl("SSP2.2050", co2_glob$variable), "SSP245 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP5.2050", co2_glob$variable), "SSP585 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP1.2050", co2_glob$variable), "SSP126 (2050)", co2_glob$Scenario ) 
co2_glob$Scenario <- ifelse(grepl("SSP3.2050", co2_glob$variable), "SSP370 (2050)", co2_glob$Scenario ) 

co2_glob <- filter(co2_glob, !is.na(Scenario))

co2_glob$q <- NA
co2_glob$q <- ifelse(grepl("q1", co2_glob$variable), "Q1", co2_glob$q ) 
co2_glob$q <- ifelse(grepl("q3", co2_glob$variable), "Q3", co2_glob$q ) 
co2_glob$q <- ifelse(is.na(co2_glob$q), "Median", co2_glob$q ) 

co2_glob$variable <- NULL

co2_glob$Region <- ifelse(co2_glob$region=="East Asia & Pacific", "CHN", NA)
co2_glob$Region <- ifelse(co2_glob$region=="Europe & Central Asia", "EU", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Latin America & Caribbean", "BRA", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Middle East & North Africa", "TUR", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="North America", "USA", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="South Asia", "IND", co2_glob$Region)
co2_glob$Region <- ifelse(co2_glob$region=="Sub-Saharan Africa", "ZAF", co2_glob$Region)

#

co2_glob <- pivot_wider(co2_glob, names_from = q, values_from = value)

#

merger <- merge(co2_glob, ci_emis, by=c("Scenario", "Region"), all.x=T)

merger <- group_by(merger,region) %>% mutate(emission_factor = ifelse(is.na(emission_factor), mean(emission_factor, na.rm=T), emission_factor)) 

merger$mtco2_q5 <- merger$Median * merger$emission_factor / 1000000000000 # convert to Mt
merger$mtco2_q1 <- merger$Q1 * merger$emission_factor / 1000000000000 # convert to Mt
merger$mtco2_q3 <- merger$Q3 * merger$emission_factor / 1000000000000 # convert to Mt

merger <- dplyr::select(merger, Scenario, ISO3, mtco2_q5, mtco2_q1, mtco2_q3)

#

merger <- merger %>% group_by(Scenario, ISO3) %>% dplyr::summarise(mtco2_q5 = mean(mtco2_q5), mtco2_q1 = mean(mtco2_q1), mtco2_q3 = mean(mtco2_q3))

mergers <- merger

merger <- merger %>%
  bind_rows(data.frame(ISO3 = "Total", Scenario=unique(merger$Scenario), mtco2_q5 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q5=sum(mtco2_q5)) %>% pull(mtco2_q5), mtco2_q1 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q1=sum(mtco2_q1)) %>% pull(mtco2_q1), mtco2_q3 = mergers %>% group_by(Scenario) %>% dplyr::summarise(mtco2_q3=sum(mtco2_q3)) %>% pull(mtco2_q3)))

#

merger$region <- NULL

colnames(merger) <- c("Scenario", "ISO3", "Mt CO2", "Mt CO2 (Q1)", "Mt CO2 (Q3)")

merger$`Mt CO2` <- round(merger$`Mt CO2`, 1)
merger$`Mt CO2 (Q1)` <- round(merger$`Mt CO2 (Q1)`, 1)
merger$`Mt CO2 (Q3)` <- round(merger$`Mt CO2 (Q3)`, 1)

merger$`Mt CO2` <- paste0(merger$`Mt CO2`, " (", merger$`Mt CO2 (Q1)`, " - ", merger$`Mt CO2 (Q3)`, ")")

merger <- dplyr::select(merger, 1,2,3)

merger <- group_by(merger, Scenario, ISO3) %>% sample_n(1)

merger <- pivot_wider(merger, names_from = Scenario, values_from = `Mt CO2`)


sink("results/graphs_tables/emissions_level_results.tex")
xtable(merger)
sink()





rm(list=setdiff(ls(), "wd")) # Removes all previously created variables
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
library(parallel)
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
require(data.table)
require(ggplot2)
require(lattice)

####

setwd(wd)

load(paste0(wd, "/results/xgboost_models_jan24.Rdata"))

##########################

library(doParallel)

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

library(pdp)
 
pdp1 <- pdp::partial(ac_model, 
        pred.var = c("mean_CDD18_db", "ln_total_exp_usd_2011"), 
        trim.outliers = TRUE, chull = FALSE, parallel = TRUE,
        grid.resolution = 30,  prob=T, paropts = list(.packages = "ranger"))

stopCluster(cluster)
gc()

pdp1$mean_CDD18_db <- exp(pdp1$mean_CDD18_db)
pdp1$ln_total_exp_usd_2011 <- exp(pdp1$ln_total_exp_usd_2011)
pdp1$yhat <- (1- pdp1$yhat) *100

####

pdp2 <- pdp::partial(ely_model, 
             pred.var = c("phat0_obs", "ln_total_exp_usd_2011", "mean_CDD18_db"), 
             trim.outliers = TRUE, chull = F, parallel = F,
             grid.resolution = 25,  paropts = list(.packages = "ranger"))

pdp2$mean_CDD18_db <- exp(pdp2$mean_CDD18_db)
pdp2$ln_total_exp_usd_2011 <- exp(pdp2$ln_total_exp_usd_2011)
pdp2$phat0_obs <- pdp2$phat0_obs*100
pdp2$yhat <- exp(pdp2$yhat)

save(pdp1, pdp2, file="pdp_outputs.Rdata")

load("pdp_outputs.Rdata")

#####

pdp1_df <- pdp1
pdp2_df <- pdp2

pdp1 <- data.table(pdp1)
pdp2 <- data.table(pdp2)

pdp2

s1_yhat_at_min_mean_CDD18_db <- pdp1[mean_CDD18_db==min(mean_CDD18_db),]$yhat
s1_yhat_at_min_ln_total_exp_usd_2011 <- pdp1[ln_total_exp_usd_2011==min(ln_total_exp_usd_2011),]$yhat

pdp1[,c('elas_cdd','elas_inc'):=
       list(0.01 * (yhat - s1_yhat_at_min_mean_CDD18_db) / (log(mean_CDD18_db) - log(min(mean_CDD18_db))),
            0.01 * (yhat - s1_yhat_at_min_ln_total_exp_usd_2011) / (log(ln_total_exp_usd_2011) - log(min(ln_total_exp_usd_2011))))]

summary(pdp1[mean_CDD18_db!=min(mean_CDD18_db) & ln_total_exp_usd_2011!=min(ln_total_exp_usd_2011),])

p_first_stage <- ggplot()+
  geom_tile(data=pdp1_df,
            aes(x=log(mean_CDD18_db),
                y=log(ln_total_exp_usd_2011),
                fill=yhat))+
  theme_classic()+
  scale_fill_gradientn(colors=c('blue','green','yellow','orange','red'), name="%")+
  theme(legend.position='bottom',legend.direction='horizontal',legend.title=element_blank())+
  xlab("Log Cooling Degree Days")+
  ylab("Log total expenditure")+
  theme(text=element_text(size=15), #change font size of all text
          axis.text=element_text(size=15), #change font size of axis text
          axis.title=element_text(size=15), #change font size of axis titles
          plot.title=element_text(size=15), #change font size of plot title
          legend.text=element_text(size=15), #change font size of legend text
          legend.title=element_text(size=15)) #change font size of legend title  

p_first_stage + ggtitle("d")

ggsave("results/graphs_tables/pdp_1.pdf", height = 5, width = 5)

###

pdp2

s2_yhat_at_min_mean_CDD18_db <- pdp2[mean_CDD18_db==min(mean_CDD18_db),]$yhat
s2_yhat_at_min_ln_total_exp_usd_2011 <- pdp2[ln_total_exp_usd_2011==min(ln_total_exp_usd_2011),]$yhat
s2_yhat_at_min_phat0_obs <- pdp2[phat0_obs==min(phat0_obs),]$yhat

pdp2[,c('elas_cdd','elas_inc','elas_ac'):=
       list((log(yhat) - log(s2_yhat_at_min_mean_CDD18_db)) / (log(mean_CDD18_db) - log(min(mean_CDD18_db))),
            (log(yhat) - log(s2_yhat_at_min_ln_total_exp_usd_2011)) / (log(ln_total_exp_usd_2011) - log(min(ln_total_exp_usd_2011))),
            (log(yhat) - log(s2_yhat_at_min_phat0_obs)) / (0.01 * (phat0_obs - min(phat0_obs))))]

summary(pdp2[mean_CDD18_db!=min(mean_CDD18_db) & 
               ln_total_exp_usd_2011!=min(ln_total_exp_usd_2011) &
               phat0_obs!=min(phat0_obs),])

pdp2_df$ln_yhat <- log(pdp2_df$yhat)
pdp2_df$ln_cdd <- log(pdp2_df$mean_CDD18_db)
pdp2_df$ln_exp <- log(pdp2_df$ln_total_exp_usd_2011)
pdp2_df$acprob <- 0.01 * pdp2_df$phat0_obs

setorder(pdp2_df,yhat)

mycols <- colorRampPalette(colors=c('blue','green','yellow','orange','red'))

# remove ticks, box around plot
# https://stat.ethz.ch/pipermail/r-help/2003-July/036826.html

pdp2_df$acprob_q <- cut(pdp2_df$acprob, (c(0, 0.25, 0.5, 0.75, 1)), labels=(c("0-25", "25-50", "50-75", "75-100")))

pdp2_df$acprob_q <- factor(pdp2_df$acprob_q, levels = rev(c("25-50", "0-25", "75-100", "50-75")))

pdf("results/graphs_tables/pdp_2.pdf", width = 8, height = 8)
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(exp(ln_yhat) ~ ln_cdd * ln_exp | acprob_q,
          data = pdp2_df, drape =T,
          xlab=list("   log Cooling Degree Days", rot=30),ylab=list("log expenditure", rot=-30) ,zlab = list("log electricity",rot=90),
          scale = list(arrows = F),col.regions=mycols(nrow(pdp2_df)),colorkey=list(space='bottom'),
          col=NA)
dev.off()

###

load("pdp_outputs.Rdata")


pdp1_income <- matrix(data=pdp1$yhat, nrow=30, ncol=30)
rownames(pdp1_income) <- unique(pdp1$mean_CDD18_db)
colnames(pdp1_income) <- unique(pdp1$ln_total_exp_usd_2011)

pdp1_income_r <- pdp1_income

for (j in 1:nrow(pdp1_income_r)){
  for (i in 2:ncol(pdp1_income_r)){
    pdp1_income_r[j,i] <-((pdp1_income[j,i]/pdp1_income[j,(i-1)])-1) / ((as.numeric(rownames(pdp1_income))[i]/as.numeric(rownames(pdp1_income))[i-1])-1) 
  }}

pdp1_income_r <- pdp1_income_r[,-1]

pdp1_income_r <- as.vector(pdp1_income_r)

###

pdp1_cdd <- matrix(data=pdp1$yhat, nrow=30, ncol=30)
rownames(pdp1_cdd) <- unique(pdp1$ln_total_exp_usd_2011)
colnames(pdp1_cdd) <- unique(pdp1$mean_CDD18_db)

#######

pdp1_cdd_r <- pdp1_cdd

for (j in 1:nrow(pdp1_cdd_r)){
  for (i in 2:ncol(pdp1_cdd_r)){
    pdp1_cdd_r[j,i] <-((pdp1_cdd[j,i]/pdp1_cdd[j,(i-1)])-1) / ((as.numeric(rownames(pdp1_cdd))[i]/as.numeric(rownames(pdp1_cdd))[i-1])-1) 
  }}

pdp1_cdd_r <- pdp1_cdd_r[,-1]

pdp1_cdd_r <- as.vector(pdp1_cdd_r)

slopes_ac <- data.frame(term=c(rep("Expenditure", 870), rep("Cooling Degree Days", 870)), estimate=c(pdp1_income_r, pdp1_cdd_r))

boxplot1 <- ggplot(slopes_ac)+
  theme_classic()+
  geom_boxplot(aes(x=term, y=estimate), fill="lightblue")+
  ylab("Elasticity of air-conditioning ownership probability")+
  xlab("Driver")+
  ggtitle("e")+
  ylim(c(-0.5, 1.5))

###

pdp2_income_elas <- pdp2

pdp2_income_elas <- arrange(pdp2_income_elas, phat0_obs, mean_CDD18_db)

pdp2_income_elas$elasticity <- NA

pdp2_income_elas <- split(pdp2_income_elas, paste0(pdp2_income_elas$phat0_obs, "_", pdp2_income_elas$mean_CDD18_db))


for(layer in 1:length(pdp2_income_elas)){
  
  for (i in 2:25){
    
    pdp2_income_elas[[layer]]$elasticity[i] <- ((pdp2_income_elas[[layer]]$yhat[i]/pdp2_income_elas[[layer]]$yhat[(i-1)])-1) / ((pdp2_income_elas[[layer]]$ln_total_exp_usd_2011[i]/pdp2_income_elas[[layer]]$ln_total_exp_usd_2011[i-1])-1) 
    
  }}

pdp2_income_elas <- bind_rows(pdp2_income_elas)

###

pdp2_cdd_elas <- pdp2

pdp2_cdd_elas <- arrange(pdp2_cdd_elas, phat0_obs, ln_total_exp_usd_2011)

pdp2_cdd_elas$elasticity <- NA

pdp2_cdd_elas <- split(pdp2_cdd_elas, paste0(pdp2_cdd_elas$phat0_obs, "_", pdp2_cdd_elas$ln_total_exp_usd_2011))


for(layer in 1:length(pdp2_cdd_elas)){
  
  for (i in 2:25){
    
    pdp2_cdd_elas[[layer]]$elasticity[i] <- ((pdp2_cdd_elas[[layer]]$yhat[i]/pdp2_cdd_elas[[layer]]$yhat[(i-1)])-1) / ((pdp2_cdd_elas[[layer]]$mean_CDD18_db[i]/pdp2_cdd_elas[[layer]]$mean_CDD18_db[i-1])-1) 
    
  }}

pdp2_cdd_elas <- bind_rows(pdp2_cdd_elas)

boxplot(pdp2_cdd_elas$elasticity)

###

pdp2_ac_elas <- pdp2

pdp2_ac_elas <- arrange(pdp2_ac_elas, ln_total_exp_usd_2011, mean_CDD18_db)

pdp2_ac_elas$elasticity <- NA

pdp2_ac_elas <- split(pdp2_ac_elas, paste0(pdp2_ac_elas$ln_total_exp_usd_2011, "_", pdp2_ac_elas$mean_CDD18_db))


for(layer in 1:length(pdp2_ac_elas)){
  
  for (i in 2:25){
    
    pdp2_ac_elas[[layer]]$elasticity[i] <- ((pdp2_ac_elas[[layer]]$yhat[i]/pdp2_ac_elas[[layer]]$yhat[(i-1)])-1) / ((pdp2_ac_elas[[layer]]$phat0_obs[i]/pdp2_ac_elas[[layer]]$phat0_obs[i-1])-1) 
    
  }}

pdp2_ac_elas <- bind_rows(pdp2_ac_elas)

############
############

slopes_ely <- data.frame(term=c(rep("Expenditure", 15625), rep("Cooling Degree Days", 15625), rep("Air-conditioning", 15625)), estimate=c(pdp2_income_elas$elasticity, pdp2_cdd_elas$elasticity, pdp2_ac_elas$elasticity))

boxplot2 <- ggplot(slopes_ely)+
  theme_classic()+
  geom_boxplot(aes(x=term, y=estimate), fill="lightblue", outlier.alpha = 0)+
  ylab("Elasticity of electricity consumption")+
  xlab("Driver")+
  ylim(c(-0.5, 1.5))

library(patchwork)

boxplot1 + boxplot2 & theme(text=element_text(size=12), #change font size of all text
                              axis.text=element_text(size=12), #change font size of axis text
                              axis.title=element_text(size=12), #change font size of axis titles
                              plot.title=element_text(size=15), #change font size of plot title
                              legend.text=element_text(size=12), #change font size of legend text
                              legend.title=element_text(size=12)) #change font size of legend title  

ggsave("results/graphs_tables/boxplots_pdp.pdf", scale=1.6, width = 7, height = 3)


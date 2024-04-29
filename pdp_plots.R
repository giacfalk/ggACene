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
                fill=0.01*yhat))+
  theme_classic()+
  scale_fill_gradientn(colors=c('blue','green','yellow','orange','red'))+
  theme(legend.position='bottom',legend.direction='horizontal',legend.title=element_blank())+
  xlab("Log CDDs")+
  ylab("Log total expenditure")

p_first_stage

ggsave("results/graphs_tables/pdp_1.png")

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

png("results/graphs_tables/pdp_2.png", width = 1200, height = 1200, res = 200)
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(as.formula(ln_yhat ~ ln_cdd + ln_exp + acprob),
          data = pdp2_df, drape =T,
          xlab="log CDDs",ylab="AC Probability   ",zlab = list("log total expenditure",rot=90),
          scale = list(arrows = F),col.regions=mycols(nrow(pdp2_df)),colorkey=list(space='bottom'),
          col=NA)
dev.off()

#####

slopes_ac <- dplyr::select(pdp1, elas_cdd, elas_inc)
slopes_ac <- reshape2::melt(slopes_ac)

boxplot1 <- ggplot(slopes_ac)+
  theme_classic()+
  geom_boxplot(aes(x=variable, y=value), fill="lightblue")+
  ylab("Elasticity of AC ownership probability")+
  xlab("Driver")+
  ggtitle("D")+
  ylim(c(-0.5, 0.5))

slopes_ely <- dplyr::select(pdp2, elas_ac, elas_cdd, elas_inc)
slopes_ely <- reshape2::melt(slopes_ely)

boxplot2 <- ggplot(slopes_ely)+
  theme_classic()+
  geom_boxplot(aes(x=variable, y=value), fill="lightblue", outlier.alpha = 0)+
  ylab("Elasticity of household electricity consumption")+
  xlab("Driver")+
  ylim(c(-0.5, 1.5))

library(patchwork)

boxplot1 + boxplot2

# ggsave("results/graphs_tables/boxplots_pdp.png", scale=1.4, width = 6, height = 3)


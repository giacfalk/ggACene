
######
# ggACene (global gridded Air Conditioning energy) projections
# Giacomo Falchetta
# 26/04/2024
# https://github.com/giacfalk/ggACene

# set working directory
wd <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/global_spline"
setwd(wd)

# prepare gridded input data
# source("survey_ac_stats_table.R")
source("prepare_gridded_input_data.R")
# source("survey_regions_map.R")

# plot and describe future evolution of drivers
setwd(wd)
source("plot-drivers.R")
source("table_drivers.R")
# source("compare_training_data_preds.R")

# train the machine learning model on the household survey data
rm(list=setdiff(ls(), "wd")); gc()
training_frac <- 1
source("rf_model.R")

# additional benchmarks for the ML model
source("pdp_plots.R")
source("shap_values.R")
source("bench_splines.R")
source("models_metrics_plots.R")

# make and validate gridded projections
source("calibrate_gridded_input_data.R")
source("spline_make_gridded_projections_rf.R")
source("validation_plot.R")

# write the gridded dataset

source("writeNCDF.R")

# generate summary figures and stats

setwd(wd)
source("summary_projections_plot.R")
source("decile_boxplots.R")
source("emissions_and_cooling_gap.R")
source("cdd_income_distribution.R")
source("extra_summary_plots.R")
source("stats_for_paper.R")


######
# ggACene (global gridded Air Conditioning energy) projections
# Giacomo Falchetta
# 29/02/2024
# https://github.com/giacfalk/ggACene

# set working directory
wd <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/global_spline"
setwd(wd)

# prepare gridded input data
source("prepare_gridded_input_data.R")

# plot and describe future evolution of drivers
setwd(wd)
source("plot-drivers.R")
source("table_drivers.R")
source("compare_training_data_preds.R")

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

# generate summary figures and stats

setwd(wd)
source("summary_projections_plot.R")
source("stats_for_paper.R")
source("cumulative_figure.R")
source("emissions_and_cooling_gap.R")

# write the gridded dataset

source("writeNCDF.R")

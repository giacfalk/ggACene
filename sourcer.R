
# set working directory
wd <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/global_spline"

setwd(wd)

source("prepare_gridded_input_data_new.R")

setwd(wd)


source("plot-drivers.R")

source("table_drivers.R")

source("compare_training_data_preds.R")

rm(list=ls(all=TRUE)); gc()
training_frac <- 1

source("rf_model.R")

# additional benchmarks

source("pdp_plots.R")
source("shap_values.R")
source("bench_splines.R")

#####

source("models_metrics_plots.R")
source("calibrate_gridded_input_data.R")
source("spline_make_gridded_projections_rf.R")
source("validation_plot.R")

#

setwd(wd)

source("summary_projections_plot.R")

source("stats_for_paper.R")

source("cumulative_figure.R")

source("emissions_and_cooling_gap.R")

source("writeNCDF.R")

#############
#############

# other child projects

# source("message_comparison/aggregate_numbers.R")

# source("combine_with_elderly_projections.R")

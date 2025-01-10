# script for "Estimating the effect of COVID-19 vaccination and prior infection
# on Ct values as a proxy of SARS-CoV-2 viral load"

# prepare -----------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")

# import ------------------------------------------------------------------

source("scripts/01_import/01_import_Ct_EVT.r")

# transform ---------------------------------------------------------------

source("scripts/02_transform/02_transform_Ct_EVT.r")

# save data ---------------------------------------------------------------

data_ct %>% saveRDS("data/data_ct_20240522.rds")

# load --------------------------------------------------------------------

data_ct <- readRDS("data/data_ct_20240522.rds")

# create output -----------------------------------------------------------

source("scripts/03_tables/03_table_one.r")

source("scripts/04_figures/04_figure_hist_N.r")

source("scripts/04_figures/04_figure_hist.r")

# models ------------------------------------------------------------------

source("scripts/05_models/05.1_univariate_lm.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_inf.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_inf_time.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_vacc.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_vacc_time.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_inf_variant.r")

source("scripts/05_models/svyglm/05.2_multivariate_IPW_svyglm_vacc_variant.r")

source("scripts/05_models/svyglm/sensitivity_analysis.R")

# save output -------------------------------------------------------------

source("scripts/05_save/05_save.r")

# other -------------------------------------------------------------------

source("scripts/05_save/05_numbers_in_text.r")

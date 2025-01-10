# variables for univariate analysis:
# - age
# - sex
#
# - symptoms (T/F, or symptoms)
# - time since symptom onset
# - testing date (month)
# - variant
#
# - vaccination status / number of vaccinations
# - time since vaccination
# - vaccine type
#
# - previous infection
# - time since previous infection
# - number of previous infections

# init --------------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")

source("scripts/05_models/functions/lr_lm_wf.R")
source("scripts/05_models/functions/tidy_ci.R")
source("scripts/05_models/functions/mod_tidy.R")
source("scripts/05_models/functions/plot_estimates_univar.R")

# import, select rows/columns  --------------------------------------------

data_ct_org <- read_rds("data/data_ct_20240522.rds") %>%
  droplevels()

vars <- c(
  "age_group10_ref",
  "sex",
  "num_lab",
  "time_since_symp_onset_cat",
  "test_month",
  "variant_short_ref"
)

# split analysis by num_lab
data_ct <- data_ct_org %>%
  select(all_of(c(c(
    "Ct ORF1ab", "Ct N", "Ct S"
  ), vars))) %>%
  rename(ct_orf1ab = `Ct ORF1ab`,
         ct_n = `Ct N`,
         ct_s = `Ct S`) %>%
  mutate(num_lab = factor(num_lab))

data_ct %>% glimpse()
data_ct %>% names()


outcomes <- data_ct %>% select(starts_with("Ct")) %>% names

# split data for testing --------------------------------------------------

# make model run faster by taking half of the data
set.seed(123)
# ct_split <- initial_split(data_ct, prop = 1/5)
#
# ct_train <- training(ct_split)
# ct_test  <- testing(ct_split)

# model pipeline ----------------------------------------------------------
univariate_wf <-
  expand_grid(outcome = outcomes, var = vars) %>%
  mutate(
    str_formula = str_c(outcome, " ~ ", var),
    model = map(
      .x = str_formula,
      .f = ~ lr_lm_wf(str_formula = .x,
                       data = data_ct)
    )
  )

univariate_tidy <- univariate_wf %>%
  mutate(
    tidy = map2(.x = model,
                .y = str_formula,
                .f = tidy_ci),
    mod_tidy = map2(.x = tidy,
                    .y = str_formula,
                    .f = mod_tidy)
  )

univariate_plot <- univariate_tidy %>%
  mutate(plot = map2(.x = mod_tidy,
                     .y = str_formula,
                     .f = plot_estimates_univar))

univar_covariates <-
  plot_grid(
    plotlist = univariate_plot %>% filter(var %in% c(
      univariate_plot$var %>% unique()
    )) %>% arrange(var) %>% pull(plot),
    ncol = 3,
    align = "v"
  )

univar_covariates %>% 
  write_rds(file = "output/plot_univar_covariates.rds")


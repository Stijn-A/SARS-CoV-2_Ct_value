source("scripts/00_prepare/00_prepare.R")

dependent_vars <- c("Ct ORF1ab", "Ct N", "Ct S", "Ct ef")
#independent_vars <- c("fct_vacc_status", "time_since_vaccination", "n_previous_infection", "time_since_previous_infection")
corrections <-
  c("age_group10",
    "sex",
    "test_month",
    "variant_short",
    "time_since_symp_onset_cat")

vars %>% names

independent_vars <-
  c(
    "fct_immune_status",
    "fct_vacc_time",
    "fct_vacc_status",
    "fct_prev_inf_time",
    "n_previous_infections",
    "age_group10",
    "sex",
    "test_month",
    "variant_short",
    "time_since_symp_onset_cat"
  )


data_ct_org <- read_rds("data/data_ct_20230816.rds") %>%
  filter(Symptoms_at_least_1,
         time_since_symp_onset <= 24) %>%
  droplevels()


data_ct %>% glimpse()
data_ct %>% names()

strata <- c("lab")

# split analysis by lab
data_ct <- data_ct_org %>%
  select(all_of(c(
    c("Ct ORF1ab", "Ct N", "Ct S", "Ct ef"), independent_vars, strata
  ))) %>%
  mutate(row_nr = row_number()) %>%
  pivot_longer(
    cols = starts_with("Ct"),
    values_to = "Ct",
    names_to = "Target",
    names_prefix = "Ct "
  ) %>%
  mutate(lab_target = paste0("Ct_", Target, ".", lab)) %>%
  filter(!is.na(Ct))  %>%
  pivot_wider(
    id_cols = c(row_nr, all_of(independent_vars)),
    names_from = lab_target,
    values_from = Ct
  ) %>%
  select(!row_nr)

outcomes <- data_ct %>% select(starts_with("Ct")) %>% names

# make model run faster by taking half of the data
set.seed(123)
# ct_split <- initial_split(data_ct, prop = 1 / 5)
# 
# ct_train <- training(ct_split)
# ct_test  <- testing(ct_split)

#
# analysis 1 vacc status x n_previous_infection, 'immune status'
#

independent_vars_1 <- c(
  "fct_immune_status",
  "age_group10",
  "sex",
  "test_month",
  "variant_short",
  "time_since_symp_onset_cat"
)

source("scripts/05_models/functions/lr_glm_wf.R")
source("scripts/05_models/functions/tidy_ci.R")
source("scripts/05_models/functions/mod_tidy.R")
source("scripts/05_models/functions/plot_estimates_multivar.R")


multivar_1_wf <-
  expand_grid(outcome = outcomes,
              var = paste0(independent_vars_1, collapse = " + ")) %>%
  mutate(
    str_formula = str_c("`", outcome, "`", " ~ ", var),
    
    model = map(
      .x = str_formula,
      .f = ~ lr_glm_wf(
        str_formula = .x,
        data = data_ct  %>% filter(!is.na(outcome) &
                                     !is.na(fct_immune_status))
      )
    )
  )

multivar_1_tidy <- multivar_1_wf %>%
  mutate(
    tidy = map2(.x = model,
                .y = str_formula,
                .f = tidy_ci),
    mod_tidy = map2(.x = tidy,
                    .y = str_formula,
                    .f = mod_tidy)
  )

multivar_1_plot <- multivar_1_tidy %>%
  mutate(plot = map2(
    .x = mod_tidy,
    .y = str_formula,
    .f = ~ plot_estimates_multivar(
      glm_tidy = .x,
      str_formula = .y,
      y_filter = "fct_immune_status"
    )
  ))

multivar_1_plot %>% saveRDS(paste0(
  "output/multivar_1_plot",
  format(now(), format = "%Y%m%d_%H%M") ,
  ".rds"
))

cowplot::plot_grid(plotlist = multivar_1_plot %>% pull(plot))

test <- multivar_1_plot %>% filter(
  str_formula ==
    "`Ct_ORF1ab.Saltro` ~ fct_immune_status + age_group10 + sex + test_month + variant_short + time_since_symp_onset_cat"
) %>%
  pull(mod_tidy)

#
# analysis 2 vacc_status x time_since_vacc AND n_previous_infection x time_since_previous_infection
#

independent_vars_vacc2 <-
  c(
    "fct_vacc_time",
    "n_previous_infections",
    "age_group10",
    "sex",
    "test_month",
    "variant_short",
    "time_since_symp_onset_cat"
  )
data_ct <- data_ct_org %>%
  select(all_of(c(
    c("Ct ORF1ab", "Ct N", "Ct S", "Ct ef"), independent_vars_vacc2, strata
  ))) %>%
  mutate(row_nr = row_number()) %>%
  pivot_longer(
    cols = starts_with("Ct"),
    values_to = "Ct",
    names_to = "Target",
    names_prefix = "Ct "
  ) %>%
  mutate(lab_target = paste0("Ct_", Target, ".", lab)) %>%
  filter(!is.na(Ct))  %>%
  pivot_wider(
    id_cols = c(row_nr, all_of(independent_vars_vacc2)),
    names_from = lab_target,
    values_from = Ct
  ) %>%
  select(!row_nr)

outcomes <- data_ct %>% select(starts_with("Ct")) %>% names

multivar_vacc2_wf <-
  expand_grid(outcome = outcomes,
              var = paste0(independent_vars_vacc2, collapse = " + ")) %>%
  #filter(outcome == "Ct_ORF1ab.Synlab") %>%
  mutate(
    str_formula = str_c("`", outcome, "`", " ~ ", var),
    
    model = map(.x = str_formula,
                .f = lr_glm_wf,
                data = data_ct  %>% filter(!is.na(outcome) &
                                             !is.na(fct_vacc_time)))
  )

multivar_vacc2_tidy <- multivar_vacc2_wf %>%
  #filter(var %in% c("age_group10", "variant_short")) %>%
  mutate(
    tidy = map2(.x = model,
                .y = str_formula,
                .f = tidy_ci),
    mod_tidy = map2(.x = tidy,
                    .y = str_formula,
                    .f = mod_tidy)
  )

multivar_vacc2_plot <- multivar_vacc2_tidy %>%
  mutate(plot = map2(
    .x = mod_tidy,
    .y = str_formula,
    .f = ~ plot_estimates_multivar(
      glm_tidy = .x,
      str_formula = .y,
      y_filter = "fct_vacc_time"
    )
  ))

cowplot::plot_grid(plotlist = multivar_vacc2_plot %>% pull(plot))
multivar_vacc2_plot %>% saveRDS(paste0(
  "output/multivar_vacc2_plot",
  format(now(), format = "%Y%m%d_%H%M") ,
  ".rds"
))


# infection time ----------------------------------------------------------


independent_vars_inf2 <-
  c(
    "fct_prev_inf_time",
    "fct_vacc_status",
    "age_group10",
    "sex",
    "test_month",
    "variant_short",
    "time_since_symp_onset_cat"
  )

data_ct <- data_ct_org %>%
  select(all_of(c(
    c("Ct ORF1ab", "Ct N", "Ct S", "Ct ef"), independent_vars_inf2, strata
  ))) %>%
  mutate(row_nr = row_number()) %>%
  pivot_longer(
    cols = starts_with("Ct"),
    values_to = "Ct",
    names_to = "Target",
    names_prefix = "Ct "
  ) %>%
  mutate(lab_target = paste0("Ct_", Target, ".", lab)) %>%
  filter(!is.na(Ct))  %>%
  pivot_wider(
    id_cols = c(row_nr, all_of(independent_vars_inf2)),
    names_from = lab_target,
    values_from = Ct
  ) %>%
  select(!row_nr)

multivar_inf2_wf <-
  expand_grid(outcome = outcomes,
              var = paste0(independent_vars_inf2, collapse = " + ")) %>%
  #filter(outcome == "Ct_ORF1ab.Synlab") %>%
  mutate(
    str_formula = str_c("`", outcome, "`", " ~ ", var),
    
    model = map(.x = str_formula,
                .f = lr_glm_wf,
                data = data_ct  %>% filter(!is.na(outcome) &
                                             !is.na(fct_prev_inf_time)))
  )

multivar_inf2_tidy <- multivar_inf2_wf %>%
  #filter(var %in% c("age_group10", "variant_short")) %>%
  mutate(
    tidy = map2(.x = model,
                .y = str_formula,
                .f = tidy_ci),
    mod_tidy = map2(.x = tidy,
                    .y = str_formula,
                    .f = mod_tidy)
  )

multivar_inf2_plot <- multivar_inf2_tidy %>%
  mutate(plot = map2(
    .x = mod_tidy,
    .y = str_formula,
    .f = ~ plot_estimates_multivar(
      glm_tidy = .x,
      str_formula = .y,
      y_filter = "fct_prev_inf_time"
    )
  ))

# mutate(
#   ylab_term = ylab_term %>% str_remove(y_filter) %>%
#     str_remove_all("`")
#   %>% factor(levels = data_ct$fct_prev_inf_time %>% levels)
# )

cowplot::plot_grid(plotlist = multivar_inf2_plot  %>%  pull(plot))

multivar_inf2_plot %>% write_rds(paste0(
  "output/multivar_inf2_plot_",
  format(now(), format = "%Y%m%d_%H%M") ,
  ".rds"
))



# analysis 3 'immune status' x time_since_last_event

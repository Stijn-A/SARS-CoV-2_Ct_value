gam_ps_vacc <-
  readRDS("cache/gam_ps_vacc.rds")

source("scripts/00_prepare/00_prepare.R")

# vacc status -------------------------------------------------------------

dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_vacc_status",
    "fct_vacc_status",
    "adjustment_prior_inf",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct_org <- read_rds("data/data_ct_20240522.rds") %>%
  filter(Symptoms_at_least_1,
         time_since_symp_onset <= 24,
         !is.na(num_vacc_time),
         # with known vaccination status and time since vaccination!is.na(num_vacc_time),
         age >= 18) %>%
  droplevels()

data_ct <-
  data_ct_org %>%
  select(all_of(c(
    independent_vars, dependent_vars, "Monsternummer"
  ))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(lab = factor(lab))

tab_ps <- as_tibble(predict(gam_ps_vacc, type = "response")) %>%
  rename(unvacc = V1,
         full_vacc = V2,
         booster = V3) %>%
  bind_cols(., data_ct %>% select(fct_vacc_status, Monsternummer)) %>%
  pivot_longer(cols = -c(fct_vacc_status, Monsternummer)) %>%
  filter(fct_vacc_status == name) %>%
  mutate(
    ps_w = 1 / value,
    # max ipws is max 500
    ps_w_c =
      if_else(ps_w > 500,
              500,
              ps_w),
    no_weight = 1
  )

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, no_weight, Monsternummer),
            by = "Monsternummer")

formula_svyglm <- formula(
  "`Ct N` ~ fct_vacc_status +
                            adjustment_prior_inf +
                            ns(age, df=10) +
                            sex +
                            ns(num_testing_date, df=10) +
                            variant_short +
                            ns(time_since_symp_onset, df=10) +
                            lab"
)
formula_svyglm_no_cov <- formula("`Ct N` ~ fct_vacc_status")

tab_mod <-
  expand_grid(
    formu = c(formula_svyglm, formula_svyglm_no_cov),
    weights = c("ps_w", "ps_w_c", "no_weight")
  ) %>%
  mutate(model = map2(.x = formu, .y = weights,
                      .f = \(.x, .y) {
                        design_ct <- svydesign(ids =  ~ 1,
                                               weights =  formula(paste0("~ ", .y)),
                                               data = data_ct)
                        logger::log_info(as.character(.x[3]))
                        logger::log_info(.y)
                        svyglm(.x,
                               design = design_ct)
                      }))

tab_result <- tab_mod %>%
  mutate(stats = map(.x = model, \(.x) {
    tidy(.x, conf.int = T)
  }),
  covars = as.character(formu) %>%
    str_remove("`Ct N` ~ "))

tofind <- paste(independent_vars, collapse = "|")

tab_sensitivity_analysis_vacc <- tab_result %>%
  select(-c(formu, model)) %>%
  unnest(stats) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = str_extract(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    level = str_remove(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    statistic = paste0(estimate, " (", conf.low, " - ", conf.high, ")"),
    .before = 3
  ) %>%
  select(`co-variables` = covars, weights, variable, level , statistic)


# vacc time -------------------------------------------------------------

independent_vars <-
  c(
    "num_vacc_time",
    "fct_vacc_time",
    "adjustment_prior_inf",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

model_num_vacc_time <-
  readRDS("cache/model_num_vacc_time.rds")

data_ct <- data_ct_org %>%
  select(all_of(c(
    independent_vars, dependent_vars, "Monsternummer"
  ))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(lab = factor(lab))


tab_ps <-
  as_tibble(predict(model_num_vacc_time, type = "response")) %>%
  rename(
    `unvacc` = V1,
    `full_vacc 7-59` = V2,
    `full_vacc 60-119` = V3,
    `full_vacc 120-179` = V4,
    `full_vacc 180-299` = V5,
    `full_vacc 300-599` = V6,
    `booster 7-59`      = V7,
    `booster 60-119`    = V8,
    `booster 120-179`   = V9,
    `booster 180-299` = V10
  ) %>%
  bind_cols(., data_ct %>% select(fct_vacc_time, Monsternummer)) %>%
  pivot_longer(cols = -c(fct_vacc_time, Monsternummer)) %>%
  filter(fct_vacc_time == name) %>%
  mutate(
    ps_w = 1 / value,
    # max ipws is max 500
    ps_w_c =
      if_else(ps_w > 500,
              500,
              ps_w),
    no_weight = 1,
    fct_vacc_time = fct_vacc_time %>%
      str_replace(pattern = " ", replacement = "\n") %>%
      factor(
        levels = c(
          "unvacc",
          "full_vacc\n7-59",
          "full_vacc\n60-119",
          "full_vacc\n120-179",
          "full_vacc\n180-299",
          "full_vacc\n300-599",
          "booster\n7-59",
          "booster\n60-119",
          "booster\n120-179",
          "booster\n180-299"
        )
      )
  )

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, no_weight, Monsternummer),
            by = "Monsternummer")

formula_svyglm_vacc_time <- formula(
  "`Ct N` ~ fct_vacc_time +
                         adjustment_prior_inf +
                         ns(age, df=10) +
                         sex +
                         ns(num_testing_date, df=10) +
                         variant_short +
                         ns(time_since_symp_onset, df=10) +
                         lab"
)
formula_svyglm_no_cov_vacc_time <- formula("`Ct N` ~ fct_vacc_time")

tab_mod <-
  expand_grid(
    formu = c(formula_svyglm_vacc_time, formula_svyglm_no_cov_vacc_time),
    weights = c("ps_w", "ps_w_c", "no_weight")
  ) %>%
  mutate(model = map2(.x = formu, .y = weights,
                      .f = \(.x, .y) {
                        design_ct <- svydesign(ids =  ~ 1,
                                               weights =  formula(paste0("~ ", .y)),
                                               data = data_ct)
                        logger::log_info(as.character(.x[3]))
                        logger::log_info(.y)
                        svyglm(.x,
                               design = design_ct)
                      }))

tab_result <- tab_mod %>%
  mutate(stats = map(.x = model, \(.x) {
    tidy(.x, conf.int = T)
  }),
  covars = as.character(formu) %>%
    str_remove("`Ct N` ~ "))

tofind <- paste(independent_vars, collapse = "|")

tab_sensitivity_analysis_vacc_time <- tab_result %>%
  select(-c(formu, model)) %>%
  unnest(stats) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = str_extract(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    level = str_remove(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    statistic = paste0(estimate, " (", conf.low, " - ", conf.high, ")"),
    .before = 3
  ) %>%
  select(`co-variables` = covars, weights, variable, level , statistic)

# previous inf ------------------------------------------------------------

independent_vars <-
  c(
    "num_previous_infection",
    "fct_previous_infection",
    "adjustment_vacc",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct <-
  data_ct_org %>% 
  select(all_of(c(
    independent_vars, dependent_vars, "Monsternummer"
  ))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(lab = factor(lab))

model_num_previous_infection <-
  readRDS(
    "cache/model_num_previous_infection.rds"
  )

tab_ps <-
  as_tibble(predict(model_num_previous_infection, type = "response")) %>%
  rename(
    `no prior infection` = V1,
    `1 prior infection` = V2,
    `2+ prior infections` = V3
  ) %>%
  bind_cols(., data_ct %>% select(fct_previous_infection, Monsternummer)) %>%
  pivot_longer(cols = -c(fct_previous_infection, Monsternummer)) %>%
  filter(fct_previous_infection == name) %>%
  mutate(
    ps_w = 1 / value,
    # max ipws is max 500
    ps_w_c =
      if_else(ps_w > 500,
              500,
              ps_w),
    no_weight = 1
  )


data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, no_weight, Monsternummer),
            by = "Monsternummer")

formula_svyglm_prev_inf <-
  formula(
    "`Ct N` ~ fct_previous_infection +
    adjustment_vacc +
                         ns(age, df=10) +
                         sex +
                         ns(num_testing_date, df=10) +
                         variant_short +
                         ns(time_since_symp_onset, df=10) +
                         lab"
  )
formula_svyglm_no_cov_prev_inf <-
  formula("`Ct N` ~ fct_previous_infection")

tab_mod <-
  expand_grid(
    formu = c(formula_svyglm_prev_inf, formula_svyglm_no_cov_prev_inf),
    weights = c("ps_w", "ps_w_c", "no_weight")
  ) %>%
  mutate(model = map2(.x = formu, .y = weights,
                      .f = \(.x, .y) {
                        design_ct <- svydesign(ids =  ~ 1,
                                               weights =  formula(paste0("~ ", .y)),
                                               data = data_ct)
                        logger::log_info(as.character(.x[3]))
                        logger::log_info(.y)
                        svyglm(.x,
                               design = design_ct)
                      }))

tab_result <- tab_mod %>%
  mutate(stats = map(.x = model, \(.x) {
    tidy(.x, conf.int = T)
  }),
  covars = as.character(formu) %>%
    str_remove("`Ct N` ~ "))

tofind <- paste(independent_vars, collapse = "|")

tab_sensitivity_analysis_inf <- tab_result %>%
  select(-c(formu, model)) %>%
  unnest(stats) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = str_extract(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    level = str_remove(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    statistic = paste0(estimate, " (", conf.low, " - ", conf.high, ")"),
    .before = 3
  ) %>%
  select(`co-variables` = covars, weights, variable, level , statistic)

# prev inf time -----------------------------------------------------------

independent_vars <-
  c(
    "num_prev_inf_time",
    "fct_prev_inf_time",
    "adjustment_vacc",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct <-
  data_ct_org %>% filter(lab %in% c("Saltro", "Synlab")) %>%
  select(all_of(c(
    independent_vars, dependent_vars, "Monsternummer"
  ))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(lab = factor(lab))

model_prev_inf_time <-
  readRDS("cache/model_prev_inf_time.rds")

tab_ps <-
  as_tibble(predict(model_prev_inf_time, type = "response")) %>%
  rename(
    `no prior infection` = V1,
    `1 prior infection\n30-89` = V2,
    `1 prior infection\n90-149` = V3,
    `1 prior infection\n150-239` = V4,
    `1 prior infection\n240-449` = V5,
    `1 prior infection\n450-849` = V6,
    `2+ prior infections\n30-239` = V7,
    `2+ prior infections\n240-749` = V8
  ) %>%
  bind_cols(., data_ct %>% select(fct_prev_inf_time, Monsternummer)) %>%
  pivot_longer(cols = -c(fct_prev_inf_time, Monsternummer)) %>%
  filter(fct_prev_inf_time == name) %>%
  mutate(
    ps_w = 1 / as.numeric(value),
    # max ipws is max 500
    ps_w_c =
      if_else(ps_w > 500,
              500,
              ps_w),
    no_weight = 1
  )

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, no_weight,  Monsternummer),
            by = "Monsternummer")

formula_svyglm_prev_inf_time <-
  formula(
    "`Ct N` ~ fct_prev_inf_time +
                                          adjustment_vacc +
                                          ns(age, df=10) +
                                          sex +
                                          ns(num_testing_date, df=10) +
                                          ns(time_since_symp_onset, df=10) +
                                          lab"
  )
formula_svyglm_no_cov_prev_inf_time <-
  formula("`Ct N` ~ fct_prev_inf_time")

tab_mod <-
  expand_grid(
    formu = c(
      formula_svyglm_prev_inf_time,
      formula_svyglm_no_cov_prev_inf_time
    ),
    weights = c("ps_w", "ps_w_c", "no_weight")
  ) %>%
  mutate(model = map2(.x = formu, .y = weights,
                      .f = \(.x, .y) {
                        design_ct <- svydesign(ids =  ~ 1,
                                               weights =  formula(paste0("~ ", .y)),
                                               data = data_ct)
                        logger::log_info(as.character(.x[3]))
                        logger::log_info(.y)
                        svyglm(.x,
                               design = design_ct)
                      }))

tab_result <- tab_mod %>%
  mutate(stats = map(.x = model, \(.x) {
    tidy(.x, conf.int = T)
  }),
  covars = as.character(formu) %>%
    str_remove("`Ct N` ~ "))

tofind <- paste(independent_vars, collapse = "|")

tab_sensitivity_analysis_inf_time <- tab_result %>%
  select(-c(formu, model)) %>%
  unnest(stats) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    variable = str_extract(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    level = str_remove(string = term, pattern = tofind) %>%
      str_remove("_short|fct_"),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    statistic = paste0(estimate, " (", conf.low, " - ", conf.high, ")"),
    .before = 3
  ) %>%
  select(`co-variables` = covars, weights, variable, level , statistic)



# combine tables ----------------------------------------------------------


table_sensitivity <- bind_rows(
  tab_sensitivity_analysis_vacc %>%
    filter(variable == "vacc_status") %>%
    mutate(
      `GLM adjustment` = if_else(
        `co-variables` == "fct_vacc_status",
        "not adjusted\nfor covariables",
        "adjusted\nfor covariables"
      ),
      `IPSW` = case_when(
        weights == "ps_w" ~ "Original weights",
        weights == "ps_w_c" ~ "Trimmed weights",
        weights == "no_weight" ~ "No weights"
      ),
      model = "Vaccination"
    ) %>%
    select(-c(variable, `co-variables`, weights)) %>%
    pivot_wider(names_from = "GLM adjustment",
                values_from = statistic),
  
  tab_sensitivity_analysis_vacc_time %>%
    filter(variable == "vacc_time") %>%
    mutate(
      `GLM adjustment` = if_else(
        `co-variables` == "fct_vacc_time",
        "not adjusted\nfor covariables",
        "adjusted\nfor covariables"
      ),
      `IPSW` = case_when(
        weights == "ps_w" ~ "Original weights",
        weights == "ps_w_c" ~ "Trimmed weights",
        weights == "no_weight" ~ "No weights"
      ),
      model = "Time since vaccination"
    ) %>%
    select(-c(variable, `co-variables`, weights)) %>%
    pivot_wider(names_from = "GLM adjustment",
                values_from = statistic),
  
  tab_sensitivity_analysis_inf %>%
    filter(variable == "previous_infection") %>%
    mutate(
      `GLM adjustment` = if_else(
        `co-variables` == "fct_previous_infection",
        "not adjusted\nfor covariables",
        "adjusted\nfor covariables"
      ),
      `IPSW` = case_when(
        weights == "ps_w" ~ "Original weights",
        weights == "ps_w_c" ~ "Trimmed weights",
        weights == "no_weight" ~ "No weights"
      ),
      model = "Prior infection"
    ) %>%
    select(-c(variable, `co-variables`, weights)) %>%
    pivot_wider(names_from = "GLM adjustment",
                values_from = statistic),
  
  
  tab_sensitivity_analysis_inf_time %>%
    filter(variable == "prev_inf_time") %>%
    mutate(
      `GLM adjustment` = if_else(
        `co-variables` == "fct_prev_inf_time",
        "not adjusted\nfor covariables",
        "adjusted\nfor covariables"
      ),
      `IPSW` = case_when(
        weights == "ps_w" ~ "Original weights",
        weights == "ps_w_c" ~ "Trimmed weights",
        weights == "no_weight" ~ "No weights"
      ),
      model = "Time since prior infection"
    ) %>%
    select(-c(variable, `co-variables`, weights)) %>%
    pivot_wider(names_from = "GLM adjustment",
                values_from = statistic)
) %>%
  select(model,
         level,
         IPSW,
         `adjusted\nfor covariables`,
         `not adjusted\nfor covariables`) %>%
  arrange(model, level)

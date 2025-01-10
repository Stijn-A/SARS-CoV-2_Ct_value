# vacc --------------------------------------------------------------------


source("scripts/00_prepare/00_prepare.R")
# import, row and column selection ----------------------------------------

dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_vacc_status",
    "adjustment_prior_inf",
    "fct_vacc_status",
    "fct_previous_infection",
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
         !is.na(variant_short),
         # with known vaccination status and time since vaccination!is.na(num_vacc_time),
         !is.na(num_vacc_time),
         age >= 18) %>%
  droplevels()

data_ct <-
  data_ct_org %>%
  select(all_of(c(
    independent_vars, dependent_vars, "Monsternummer"
  ))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(lab = factor(lab))

logger::log_info("Nr of records: ", data_ct %>% nrow)

gam_ps_vacc <- read_rds("cache/gam_ps_vacc.rds")

tab_ps <- as_tibble(predict(gam_ps_vacc, type = "response")) %>%
  rename(unvacc = V1,
         full_vacc = V2,
         booster = V3) %>%
  bind_cols(., data_ct %>% select(fct_vacc_status, Monsternummer)) %>%
  pivot_longer(cols = -c(fct_vacc_status, Monsternummer)) %>%
  filter(fct_vacc_status == name) %>%
  mutate(ps_w = 1 / value,
         # max ipws is max 500
         ps_w_c =
           if_else(ps_w > 500,
                   500,
                   ps_w))

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
            by = "Monsternummer")

max(tab_ps$ps_w) / sum(tab_ps$ps_w) * 100

# wGLM survey -------------------------------------------------------------


design_ct <- svydesign(ids =  ~ 1,
                       weights =  ~ ps_w_c,
                       data = data_ct)


models_targets <- c("Ct ORF1ab") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_vacc_status +
                adjustment_prior_inf +
                ns(age, df=10) +
                sex +
                ns(num_testing_date, df=10) +
                ns(time_since_symp_onset, df=10) +
                lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_ct)
  })

# plot result -------------------------------------------------------------

plot_outcomes_fct_vacc_status <- models_targets %>%
  map(.f = \(x) {
    form <-
      x %>% tidytable::pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_vacc_status")) %>%
      add_row(term = "unvacc") %>%
      mutate(
        term = str_remove(term, "fct_vacc_status") %>%
          factor(levels = data_ct$fct_vacc_status %>% levels),
        label = if_else(term == "unvacc", "reference", NA),
        term_c = fct_recode(
          term,
          "primary\nvaccination" = "full_vacc",
          "booster\nvaccination" = "booster",
          "unvaccinated" = "unvacc"
        ) %>%
          factor(
            .,
            levels = c(
              "booster\nvaccination",
              "primary\nvaccination",
              "unvaccinated"
            )
          )
      ) %>%
      ggplot(aes(
        estimate,
        term_c,
        xmin = conf.low,
        xmax = conf.high,
        height = 0,
        label = label
      )) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 4) +
      geom_errorbarh() +
      geom_text(x = 1.4) +
      scale_x_continuous(breaks = seq(-1, 3, 1)) +
      coord_cartesian(xlim = c(-1, 3)) +
      labs(y = "", x = "") +
      theme(legend.position = "none", text = element_text(size = 20))
  }) %>%
  cowplot::plot_grid(plotlist = .)

# save --------------------------------------------------------------------

plot_outcomes_fct_vacc_status %>%
  write_rds(file = "output/plot_outcomes_fct_vacc_status_orf1ab.rds")

# vacc time ---------------------------------------------------------------


source("scripts/00_prepare/00_prepare.R")

# import, row and column selection ----------------------------------------

dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_vacc_time",
    "fct_vacc_time",
    "fct_previous_infection",
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
         !is.na(variant_short),
         #test_month %in% seq(as_date("2021-03-01"), as_date("2022-07-01"), by = "month"),
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

logger::log_info("Nr of records: ", data_ct %>% nrow)

model_num_vacc_time <- read_rds("cache/model_num_vacc_time.rds")


fct_levels <- data_ct$fct_vacc_time %>% levels

tab_ps <-
  as_tibble(predict(model_num_vacc_time, type = "response")) %>%
  rename_at(
    vars(1:length(fct_levels)), 
    ~fct_levels
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
    fct_vacc_time = fct_vacc_time %>%
      str_replace(pattern = " ", replacement = "\n") %>%
      factor(
        levels = fct_levels
      )
  )

rm(fct_levels)

#tab_ps$ps_w %>% sum /
sum(tab_ps$ps_w) / 100

max(tab_ps$ps_w_c) / sum(tab_ps$ps_w_c) * 100

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
            by = "Monsternummer")

# wGLM survey -------------------------------------------------------------

design_ct <- svydesign(ids =  ~ 1,
                       weights =  ~ ps_w_c,
                       data = data_ct)

models_targets <- c("Ct ORF1ab") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_vacc_time +
                adjustment_prior_inf +
                ns(age, df=10) +
                sex +
                ns(num_testing_date, df=10) +
                ns(time_since_symp_onset, df=10) +
                lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_ct)
  })

# plot result -------------------------------------------------------------

plot_outcomes_fct_vacc_time <- models_targets %>%
  map(.f = \(x) {
    form <-
      x %>% tidytable::pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_vacc_time")) %>%
      add_row(term = "unvacc") %>%
      mutate(
        term = str_remove(term, "fct_vacc_time") %>%
          str_replace(pattern = " ", replacement = "\n") %>%
          factor(
            levels = rev(c(
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
            ))
          ),
        term_c = fct_recode(
          term,
          "unvaccinated" = "unvacc",
          "primary\n7-59" = "full_vacc\n7-59",
          "primary\n60-119" = "full_vacc\n60-119",
          "primary\n120-179" = "full_vacc\n120-179",
          "primary\n180-299" = "full_vacc\n180-299",
          "primary\n300-599" = "full_vacc\n300-599",
          "booster\n7-59" = "booster\n7-59",
          "booster\n60-119" = "booster\n60-119",
          "booster\n120-179" = "booster\n120-179",
          "booster\n180-299" = "booster\n180-299"
        ),
        label = if_else(term == "unvacc", "reference", NA)
      ) %>%
      ggplot(aes(
        estimate,
        term_c,
        xmin = conf.low,
        xmax = conf.high,
        height = 0,
        label = label
      )) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 4) +
      geom_errorbarh() +
      geom_text(x = 1.4) +
      scale_x_continuous(breaks = seq(-1, 3, 1)) +
      coord_cartesian(xlim = c(-1, 3)) +
      labs(y = "") +
      theme(
        legend.position = "none",
        text = element_text(size = 20),
        axis.title.x = element_blank()
      )
    # ggtitle(form)
  }) %>%
  cowplot::plot_grid(plotlist = .)

# save --------------------------------------------------------------------
plot_outcomes_fct_vacc_time %>% write_rds("output/plot_outcomes_fct_vacc_time_orf1ab.rds")


# prior inf ---------------------------------------------------------------

# init  -------------------------------------------------------------------
source("scripts/00_prepare/00_prepare.R")

# import, row and column selection ----------------------------------------
dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_previous_infection",
    "adjustment_vacc",
    "fct_previous_infection",
    "fct_vacc_status",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct_org <- read_rds("data/data_ct_20240522.rds") %>%
  filter(
    Symptoms_at_least_1,
    time_since_symp_onset <= 24,
    !is.na(variant_short),
    # with known vaccination status and time since vaccination
    !is.na(num_vacc_time),
    age >= 18
  ) %>%
  droplevels()

data_ct <-
  data_ct_org %>% 
  select(all_of(c(independent_vars, dependent_vars, "Monsternummer"))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct N`)) %>%
  mutate(
    lab = factor(lab)
  )
logger::log_info("Nr of records: ", data_ct %>% nrow)
# propensity score weighting ----------------------------------------------
#https://cran.r-project.org/web/packages/twang/vignettes/mnps.pdf

model_num_previous_infection <- read_rds("cache/model_num_previous_infection.rds")

tab_ps <- as_tibble(predict(model_num_previous_infection, type = "response")) %>% 
  rename(`no prior infection` = V1, 
         `1 prior infection` = V2, 
         `2+ prior infections` = V3) %>% 
  bind_cols(., data_ct %>% select(fct_previous_infection, Monsternummer)) %>% 
  pivot_longer(cols = -c(fct_previous_infection, Monsternummer)) %>% 
  filter(fct_previous_infection == name) %>% 
  mutate(ps_w = 1 / value,
         # max ipws is max 500
         ps_w_c =  
           if_else(
             ps_w > 500,
             500,
             ps_w
           ))

data_ct <- data_ct %>% 
  left_join(
    tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
    by = "Monsternummer"
  )

# wGLM survey -------------------------------------------------------------

design_ct <- svydesign(ids =  ~ 1,
                       weights =  ~ ps_w_c,
                       data = data_ct)


models_targets <- c("Ct N") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_previous_infection +
            adjustment_vacc +
            ns(age, df=10) +
            sex +
            ns(num_testing_date, df=10) +
            ns(time_since_symp_onset, df=10) +
            lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_ct)
  })

# plot result -------------------------------------------------------------

plot_outcomes_fct_previous_infection <- models_targets %>%
  map(.f = \(x) {
    form <- x %>% tidytable::pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_previous_infection")) %>%
      add_row(term = "no prior infection") %>%
      mutate(
        term = str_remove(term, "fct_previous_infection") %>%
          factor(levels = rev(data_ct$fct_previous_infection %>% levels)),
        label = if_else(term == "no prior infection", "reference", NA)
      ) %>%
      ggplot(aes(
        estimate,
        term,
        xmin = conf.low,
        xmax = conf.high,
        height = 0,
        label = label
      )) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 4) +
      geom_errorbarh() +
      geom_text(x = 1.4) +
      scale_x_continuous(breaks = seq(-1,3,1)) +
      coord_cartesian(xlim = c(-1,3)) +
      labs( y = "", x = "") +
      theme(legend.position = "none", text = element_text(size = 20))
  }) %>%
  cowplot::plot_grid(plotlist = .)


# save --------------------------------------------------------------------

plot_outcomes_fct_previous_infection %>% 
  write_rds(file = "output/plot_outcomes_fct_previous_infection_orf1ab.rds")


# time prior inf ----------------------------------------------------------


source("scripts/00_prepare/00_prepare.R")
library(survey)
library(splines)

# import, row and column selection ----------------------------------------

dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_prev_inf_time",
    "adjustment_vacc",
    "fct_prev_inf_time",
    "fct_vacc_status",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct_org <- read_rds("data/data_ct_20240522.rds") %>%
  filter(
    Symptoms_at_least_1,
    time_since_symp_onset <= 24,
    # with known vaccination status and time since vaccination
    !is.na(num_vacc_time),
    age >= 18
  ) %>%
  droplevels()

data_ct <-
  data_ct_org %>% 
  select(all_of(c(independent_vars, dependent_vars, "Monsternummer"))) %>%
  drop_na(!c(`Ct ORF1ab`,  `Ct N`)) %>%
  mutate(
    lab = factor(lab)
  )

logger::log_info("Nr of records: ", data_ct %>% nrow)

model_prev_inf_time <- read_rds("cache/model_prev_inf_time.rds")


# propensity score weighting ----------------------------------------------
#https://cran.r-project.org/web/packages/twang/vignettes/mnps.pdf

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
  mutate(ps_w = 1 / as.numeric(value),
         # max ipws is max 500
         ps_w_c =  
           if_else(
             ps_w > 500,
             500,
             ps_w
           ))

# boxplot weights ---------------------------------------------------------
data_ct <- data_ct %>% 
  left_join(
    tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
    by = "Monsternummer"
  )

# wGLM survey -------------------------------------------------------------

design_ct <- svydesign(ids =  ~ 1,
                       weights =  ~ ps_w_c,
                       data = data_ct)

models_targets <- c("Ct N") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_prev_inf_time +
            adjustment_vacc +
            ns(age, df=10) +
            sex +
            ns(num_testing_date, df=10) +
            ns(time_since_symp_onset, df=10) +
            lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_ct)
  })

# plot result -------------------------------------------------------------

plot_outcomes_fct_prev_inf_time <- models_targets %>%
  map(.f = \(x) {
    form <- x %>% tidytable::pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_prev_inf_time")) %>%
      add_row(term = "no prior infection") %>%
      mutate(
        term = str_remove(term, "fct_prev_inf_time") %>%
          factor(levels = rev(data_ct$fct_prev_inf_time %>% levels)),
        label = if_else(term == "no prior infection", "reference", NA)
      ) %>%
      ggplot(aes(
        estimate,
        term,
        xmin = conf.low,
        xmax = conf.high,
        height = 0,
        label = label
      )) +
      geom_point() +
      geom_vline(xintercept = 0, lty = 4) +
      geom_errorbarh() +
      geom_text(x = 1.4) +
      scale_x_continuous(breaks = seq(-1,3,1)) +
      coord_cartesian(xlim = c(-1,3)) +
      labs( y = "") +
      theme(legend.position = "none", text = element_text(size = 20),
            axis.title.x = element_blank())
  }) %>%
  cowplot::plot_grid(plotlist = .)

# save --------------------------------------------------------------------
plot_outcomes_fct_prev_inf_time %>% 
  write_rds(file = "output/plot_outcomes_fct_prev_inf_time_orf1ab.rds")


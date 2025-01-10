model_num_previous_infection <- readRDS("~/01_projects/02_sars_cov_2/ct_values_evt/cache/model_num_previous_infection.rds")

# init  -------------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")
library(survey)
library(splines)
# import, row and column selection ----------------------------------------
#
# log schaal?
#
#

dependent_vars <- c("Ct ORF1ab", "Ct N")
independent_vars <-
  c(
    "num_previous_infection",
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

max(tab_ps$ps_w) / sum(tab_ps$ps_w) * 100

# wGLM survey -------------------------------------------------------------

data_ct <- data_ct %>% 
  mutate(
    omicron = if_else(variant_short %in% c("Pre-VOC", 
                                           "Alpha",
                                           "Gamma",
                                           "Delta"),
                      "Non-Omicron",
                      "Omicron")
  )

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
                         fct_vacc_status +
                         omicron +
                         fct_previous_infection:omicron +
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

models_targets[[1]] %>% 
  tidy(conf.int = T) %>% 
  #select(term, estimate, conf.low, conf.high) %>% 
  # mutate(
  #   across(where(is.numeric), ~ round(., 2))) %>% 
  view()


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
          factor(levels = data_ct$fct_previous_infection %>% levels),
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

models_targets[[1]] %>% 
  tidy(conf.int = T) %>% 
  select(term, estimate, conf.low, conf.high, p.value) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2))) %>% 
  saveRDS("cache/mod_outcome_inf_variant_interaction.rds")


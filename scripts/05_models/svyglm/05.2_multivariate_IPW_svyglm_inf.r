# init  -------------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")

# import, row and column selection ----------------------------------------
#
# log schaal?
#
#

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

continuous_spline <- "ps"
ctrl <- gam.control(trace = TRUE)

logger::log_info("start fitting ps...")
model_num_previous_infection <- gam(
  list(
    num_previous_infection ~ 
      adjustment_vacc +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab,
    ~ adjustment_vacc +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab
  ),
  family = mgcv::multinom(K = 2),
  data = data_ct,
  method = "REML",
  control = ctrl,
)

gam.check(model_num_previous_infection)
summary(model_num_previous_infection)

logger::log_info("end fitting ps...")

dir.create("cache", showWarnings = F)
model_num_previous_infection %>% write_rds("cache/model_num_previous_infection.rds")

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

# boxplot weights ---------------------------------------------------------

boxplot_fct_previous_infection_large <- tab_ps %>%
  ggplot(aes(y = ps_w, x = fct_previous_infection)) +
  geom_boxplot(width = 0.5) +
  stat_summary(fun = median, fun.ymax = length,
               geom = "text", aes(label = ..ymax..), vjust = -1.8) +
  geom_hline(yintercept = 500, linetype = 2, color = "red") +
  coord_flip(ylim =  c(-1, 1000)) +
  theme_minimal()

boxplot_fct_previous_infection_small <- 
  boxplot_fct_previous_infection_large + 
  coord_flip(ylim =  c(-1, 20))

boxplot_fct_previous_infection <- cowplot::plot_grid(
  boxplot_fct_previous_infection_large, 
  boxplot_fct_previous_infection_small,
  nrow = 2
)

data_ct <- data_ct %>% 
  left_join(
    tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
    by = "Monsternummer"
  )

max(tab_ps$ps_w) / sum(tab_ps$ps_w) * 100

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

models_targets[[1]] %>% 
  tidy(conf.int = T) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2))) %>% 
  saveRDS("cache/mod_outcome_inf_ps.rds")


boxplot_fct_previous_infection %>%
  write_rds(file = "output/boxplot_fct_previous_infection_ps.rds")

plot_outcomes_fct_previous_infection %>% 
  write_rds(file = "output/plot_outcomes_fct_previous_infection_N.rds")

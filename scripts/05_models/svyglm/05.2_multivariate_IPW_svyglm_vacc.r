# init  -------------------------------------------------------------------

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
# propensity score weighting ----------------------------------------------

continuous_spline <- "ps"

ctrl <- gam.control(trace = TRUE)

logger::log_info("start fitting ps...")

gam_ps_vacc <- gam(
  list(
    num_vacc_status ~
      adjustment_prior_inf +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab,
    ~ adjustment_prior_inf +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab
  ),
  family = mgcv::multinom(K = 2),
  data = data_ct,
  method = "REML",
  control = ctrl
)
logger::log_info("end fitting ps...")

dir.create("cache", showWarnings = F)
gam_ps_vacc %>% write_rds("cache/gam_ps_vacc.rds")

gam.check(gam_ps_vacc)
summary(gam_ps_vacc)

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

# boxplot weights ---------------------------------------------------------


boxplot_fct_vacc_status_large <- tab_ps %>%
  ggplot(aes(y = ps_w, x = fct_vacc_status)) +
  geom_boxplot(width = 0.5) +
  stat_summary(
    fun = median,
    fun.max = length,
    geom = "text",
    aes(label = ..ymax..),
    vjust = -2
  ) +
  geom_hline(yintercept = 500,
             linetype = 2,
             color = "red") +
  coord_flip(ylim =  c(-1, 1000)) +
  theme_minimal()

boxplot_fct_vacc_status_small <-
  boxplot_fct_vacc_status_large +
  coord_flip(ylim =  c(-1, 20))

boxplot_fct_vacc_status <- cowplot::plot_grid(boxplot_fct_vacc_status_large,
                                              boxplot_fct_vacc_status_small,
                                              nrow = 2)

data_ct <- data_ct %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
            by = "Monsternummer")

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

models_targets[[1]] %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>% saveRDS("cache/mod_outcome_vacc.rds")

boxplot_fct_vacc_status %>%
  write_rds(file = "output/boxplot_fct_vacc_status.rds")

plot_outcomes_fct_vacc_status %>%
  write_rds(file = "output/plot_outcomes_fct_vacc_status_N.rds")

# niet log getransformeerd
# treatment levels reduceren

# gam multinom schatten
# list met aantal nomials, outcome 1x, functie = K
# splines voor continue vars, random effect voor de rest?

# gewichten vergelijken met boxplot

# svy glm met ns() voor outcome ~ treatment + cov



# init  -------------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")
library(twang)
library(survey)
library(splines)
# import, row and column selection ----------------------------------------
#
# log schaal?
#
#

dependent_vars <- c("Ct ORF1ab", "Ct N", "Ct S", "Ct ef")
independent_vars <-
  c(
    "num_vacc_status",
    "fct_vacc_status",
    "fct_previous_infection",
    "age",
    "sex",
    "num_testing_date",
    "variant_short",
    "time_since_symp_onset",
    "lab"
  )

data_ct_org <- read_rds("data/data_ct_20230907.rds") %>%
  filter(
    Symptoms_at_least_1,
    time_since_symp_onset <= 24,!is.na(fct_vacc_status),
    test_month %in% seq(as_date("2021-03-01"), as_date("2022-07-01"), by = "month"),
    age >= 18
  ) %>%
  droplevels()

data_ct <-
  data_ct_org %>% filter(lab %in% c("Saltro", "Synlab")) %>%
  select(all_of(c(independent_vars, dependent_vars))) %>%
  drop_na(!c(`Ct ORF1ab`, `Ct S`, `Ct N`, `Ct ef`)) %>%
  mutate(
    lab = factor(lab)
  )

# propensity score weighting ----------------------------------------------

continuous_spline <- "cr"
factor_spline <- "fs"
ctrl <- gam.control(trace = TRUE)

logger::log_info("start fitting ps...")


gam_ps_full_vacc <- gam(
    num_vacc_status ~ 
      fct_previous_infection +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(variant_short, bs = factor_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab,
  family = binomial,
  data = data_ct %>% filter(num_vacc_status %in% c(0, 1)),
  method = "REML",
  optimizer = c("outer", "bfgs"),
  control=ctrl
)

gam_ps_booster <- gam(
  num_vacc_status ~ 
    fct_previous_infection +
    s(age, bs = continuous_spline) +
    sex +
    s(num_testing_date, bs = continuous_spline) +
    s(variant_short, bs = factor_spline) +
    s(time_since_symp_onset, bs = continuous_spline) +
    lab,
  family = binomial,
  data = data_ct %>% filter(num_vacc_status %in% c(0, 2)),
  method = "REML",
  optimizer = c("outer", "bfgs"),
  control=ctrl
)

logger::log_info("end fitting ps...")

dir.create("cache", showWarnings = F)
gam_ps_vacc %>% write_rds("cache/gam_ps_vacc.rds")


tab_ps_full_vacc <- as_tibble(predict(gam_ps_full_vacc, type = "response")) %>% 
  rename(ps = value) %>% 
  bind_cols(.,
            data_ct %>% filter(num_vacc_status %in% c(0, 1))) %>%
  mutate(
    ipsw = if_else(fct_vacc_status == "full_vacc",
                  1 / ps, 
                  1 / (1 - ps)), # if unvacc
    comparison = "unvacc vs full vacc"
    )

tab_ps_booster <- as_tibble(predict(gam_ps_booster, type = "response")) %>% 
  rename(ps = value) %>% 
  bind_cols(.,
            data_ct %>% filter(num_vacc_status %in% c(0, 2))) %>%
  mutate(
    ipsw = if_else(fct_vacc_status == "booster",
                  1 / ps, 
                  1 / (1 - ps)), # if unvacc
    comparison = "unvacc vs booster"
    )


boxplot_fct_vacc_status_large <- 
  bind_rows(tab_ps_full_vacc, tab_ps_booster) %>%
  ggplot(aes(y = ipsw, x = fct_vacc_status)) +
  geom_boxplot(width = 0.5) +
  stat_summary(fun = median, fun.max = length,
               geom = "text", aes(label = ..ymax..), vjust = -2) +
  coord_flip(ylim =  c(-1, 1000)) +
  facet_wrap(vars(comparison), scales = "free_y")

boxplot_fct_vacc_status_small <- 
  boxplot_fct_vacc_status_large + 
  coord_flip(ylim =  c(-1, 20))

boxplot_fct_vacc_status <- cowplot::plot_grid(
  boxplot_fct_vacc_status_large, 
  boxplot_fct_vacc_status_small,
  nrow = 2
)

data_ct <- data_ct %>% 
  bind_cols(
    tab_ps_booster %>% select(ipsw)
  )

max(tab_ps$ps_w) / sum(tab_ps$ps_w) * 100

design_full_vacc <- svydesign(ids =  ~ 1,
                       weights =  ~ ipsw,
                       data = tab_ps_full_vacc)

design_booster <- svydesign(ids =  ~ 1,
                              weights =  ~ ipsw,
                              data = tab_ps_booster)


models_full_vacc <- c("Ct N", "Ct ORF1ab", "Ct S") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_vacc_status +
                         fct_previous_infection +
                         ns(age) +
                         sex +
                         ns(num_testing_date) +
                         variant_short +
                         ns(time_since_symp_onset) +
                         lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_full_vacc)
  })

models_booster <- c("Ct N", "Ct ORF1ab", "Ct S") %>%
  map(.f = \(x) {
    formula <- formula(
      paste0(
        "`",
        x,
        "`",
        " ~ fct_vacc_status +
                         fct_previous_infection +
                         ns(age) +
                         sex +
                         ns(num_testing_date) +
                         variant_short +
                         ns(time_since_symp_onset) +
                         lab"
      )
    )
    print(formula)
    svyglm(formula,
           design = design_booster)
  })

plot_outcomes_full_vacc <- models_full_vacc %>%
  map(.f = \(x) {
    form <- x %>% pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_vacc_status")) %>%
      add_row(term = "unvacc") %>%
      mutate(
        term = str_remove(term, "fct_vacc_status") %>%
          factor(levels = data_ct$fct_vacc_status %>% levels),
        label = if_else(term == "unvacc", "reference", NA)
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
      theme(legend.position = "none") +
      ggtitle(form)
  }) %>%
  cowplot::plot_grid(plotlist = .)

plot_outcomes_booster <- models_booster %>%
  map(.f = \(x) {
    form <- x %>% pull(formula) %>% as.character() %>% pluck(2) %>%
      str_remove_all("`")
    x %>% tidy(conf.int = T) %>%
      filter(str_detect(term, "fct_vacc_status")) %>%
      add_row(term = "unvacc") %>%
      mutate(
        term = str_remove(term, "fct_vacc_status") %>%
          factor(levels = data_ct$fct_vacc_status %>% levels),
        label = if_else(term == "unvacc", "reference", NA)
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
      theme(legend.position = "none") +
      ggtitle(form)
  }) %>%
  cowplot::plot_grid(plotlist = .)



boxplot_fct_vacc_status %>%
  write_rds(file = "output/boxplot_fct_vacc_status.rds")

plot_outcomes_fct_vacc_status %>% 
  write_rds(file = "output/plot_outcomes_fct_vacc_status.rds")

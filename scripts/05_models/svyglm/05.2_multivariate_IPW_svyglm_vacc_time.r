# niet log getransformeerd
# aantal treatment levels reduceren

# gam multinom schatten
# list met aantal nomials, outcome 1x, aantal formulas = K
# splines voor continue vars, random effect voor de rest?

# gewichten vergelijken met boxplot

# svy glm met ns() voor outcome ~ treatment + cov

# init  -------------------------------------------------------------------

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

# propensity score weighting ----------------------------------------------
#https://cran.r-project.org/web/packages/twang/vignettes/mnps.pdf
continuous_spline <- "cr"
ctrl <- gam.control(trace = TRUE)

logger::log_info("start fitting ps...")
model_num_vacc_time <- gam(
  list(
    num_vacc_time ~
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
      lab,
    ~ adjustment_prior_inf +
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
      lab,
    ~ adjustment_prior_inf +
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
      lab,
    ~ adjustment_prior_inf +
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
      lab,
    ~ adjustment_prior_inf +
      s(age, bs = continuous_spline) +
      sex +
      s(num_testing_date, bs = continuous_spline) +
      s(time_since_symp_onset, bs = continuous_spline) +
      lab
  ),
  family = mgcv::multinom(K = 9),
  data = data_ct,
  optimizer = c("outer", "bfgs"),
  method = "REML",
  control = ctrl
)

logger::log_info("end fitting ps...")

dir.create("cache", showWarnings = F)

model_num_vacc_time %>% write_rds("cache/model_num_vacc_time.rds")

gam.check(model_num_vacc_time)
summary(model_num_vacc_time)

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

boxplot_vacc_time_large <- tab_ps %>%
  ggplot(aes(y = ps_w, x = fct_vacc_time)) +
  geom_boxplot(width = 0.5) +
  stat_summary(
    fun = median,
    fun.ymax = length,
    geom = "text",
    aes(label = ..ymax..),
    vjust = -0.4
  ) +
  geom_hline(yintercept = 500,
             linetype = 2,
             color = "red") +
  coord_flip(ylim =  c(-1, 1000)) +
  theme_minimal()

boxplot_vacc_time_small <-
  boxplot_vacc_time_large +
  coord_flip(ylim =  c(-1, 20))

boxplot_vacc_time <- cowplot::plot_grid(boxplot_vacc_time_large,
                                        boxplot_vacc_time_small,
                                        nrow = 2)
data_ct <- data_ct %>%
  #select(-fct_vacc_time) %>%
  left_join(tab_ps %>% select(ps_w, ps_w_c, Monsternummer),
            by = "Monsternummer")


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
models_targets[[1]] %>%
  tidy(conf.int = T) %>%
  select(term, estimate, conf.low, conf.high) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  saveRDS("cache/mod_outcome_vacc_time.rds")

boxplot_vacc_time %>% write_rds("output/boxplot_vacc_time.rds")

plot_outcomes_fct_vacc_time %>% write_rds("output/plot_outcomes_fct_vacc_time_N.rds")

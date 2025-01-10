# niet log getransformeerd
# treatment levels reduceren

# gam multinom schatten
# list met aantal nomials, outcome 1x, functie = K
# splines voor continue vars, random effect voor de rest?

# gewichten vergelijken met boxplot

# svy glm met ns() voor outcome ~ treatment + cov

# init  -------------------------------------------------------------------

source("scripts/00_prepare/00_prepare.R")
library(survey)
library(splines)

# import, row and column selection ----------------------------------------

dependent_vars <- c("Ct ORF1ab", "Ct N")
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

logger::log_info("start fitting ps...")

fun_gam <- function(data, variant) {
  continuous_spline <- "ps"
  ctrl <- gam.control(trace = TRUE)
  logger::log_info("Nr of records: ", data %>% nrow)
  if (variant == "Alpha") {
  gam_ps_vacc_variant <- gam(
    list(
      num_previous_infection ~ 
        adjustment_vacc +
        s(age, bs = continuous_spline) +
        sex +
        s(num_testing_date, bs = continuous_spline) +
        s(time_since_symp_onset, bs = continuous_spline)
    ),
    family = mgcv::multinom(K = 1),
    data = data,
    method = "REML",
    control=ctrl
  )
  } else if (variant == "Omicron BA.4/5") {
    gam_ps_vacc_variant <- gam(
      list(
        num_previous_infection ~ 
          adjustment_vacc +
          s(age, bs = continuous_spline) +
          sex +
          s(num_testing_date, bs = continuous_spline) +
          s(time_since_symp_onset, bs = continuous_spline)
        ,
        ~ adjustment_vacc +
          s(age, bs = continuous_spline) +
          sex +
          s(num_testing_date, bs = continuous_spline) +
          s(time_since_symp_onset, bs = continuous_spline) 
      ),
      family = mgcv::multinom(K = 2),
      data = data,
      method = "REML",
      control=ctrl
    )
  } else {
    gam_ps_vacc_variant <- gam(
      list(
        num_previous_infection ~ 
          adjustment_vacc +
          s(age, bs = continuous_spline) +
          sex +
          s(num_testing_date, bs = continuous_spline) +
          s(time_since_symp_onset, bs = continuous_spline) +
          lab
        ,
        ~ adjustment_vacc +
          s(age, bs = continuous_spline) +
          sex +
          s(num_testing_date, bs = continuous_spline) +
          s(time_since_symp_onset, bs = continuous_spline) +
          lab
      ),
      family = mgcv::multinom(K = 2),
      data = data,
      method = "REML",
      control=ctrl
    )
  }
  
  logger::log_info("end fitting ps...")
  
  return(gam_ps_vacc_variant)
}

fun_ips <- function(gam, data) {
  col_names <- c(
    `no prior infection` = "V1", 
    `1 prior infection` = "V2", 
    `2+ prior infections` = "V3"
  )
  tab_ps <- as_tibble(predict(gam, type = "response")) %>%
    rename(any_of(col_names)) %>%
    bind_cols(.,
              data %>% 
                select(fct_previous_infection)) %>%
    pivot_longer(cols = !fct_previous_infection) %>%
    filter(fct_previous_infection == name) %>%
    mutate(ps_w = 1 / value,
           # max ipws is max 500
           ps_w_c =  
             if_else(
               ps_w > 500,
               500,
               ps_w
             ))
  
}

fun_svyglm <- function(data_ips, variant) {
  design_ct <- svydesign(ids =  ~ 1,
                         weights =  ~ ps_w_c,
                         data = data_ips)
  
  
  models_targets <- c("Ct N") %>%
    map(.f = \(x) {
      if (variant %in% c("Alpha", "Omicron BA.4/5")) {
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
                         ns(time_since_symp_onset, df=10)"
          )
        )
      } else {
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
            # lab weg
          )
        )
      }
      print(formula)
      svyglm(formula,
             design = design_ct)
    })
}

tab_inf_variant  <-
  tibble(variant = c(
    "Alpha", 
    "Delta",
    "Omicron BA.1",
    "Omicron BA.2",
    "Omicron BA.4/5"
  )) %>%
  mutate(
    data = map(.x = variant, 
               .f = ~data_ct %>% filter(variant_short == .x)),
    gam = map2(.x = data,
               .y = variant,
              .f = ~ fun_gam(data = .x, variant = .y)))

tab_inf_variant  %>% write_rds("cache/tab_inf_variant.rds")

tab_inf_variant  <- tab_inf_variant  %>%
  mutate(
    ips = map2(
      .x = gam,
      .y = data,
      .f = ~ fun_ips(gam = .x, data = .y)
    ),
    data_ips = map2(
      .x = data,
      .y = ips,
      .f = ~ .x %>%
        bind_cols(.y %>% select(ps_w, ps_w_c))
    )
  )

# boxplot weights ---------------------------------------------------------

for (i in c("Alpha", "Delta", "Omicron BA.1", "Omicron BA.2", "Omicron BA.4/5")) {
  boxplot_fct_inf_status <- tab_inf_variant  %>%
    filter(variant == i) %>%
    unnest(data_ips) %>%
    ggplot(aes(y = ps_w, x = fct_previous_infection)) +
    geom_boxplot(width = 0.5) +
    geom_hline(yintercept = 500, linetype = 2, color = "red") +
    stat_summary(
      fun = median,
      fun.max = length,
      geom = "text",
      aes(label = ..ymax..),
      vjust = -0.4
    ) +
    coord_flip(ylim =  c(-1, 20)) +
    ggtitle(i)
  
  assign(paste0("boxplot_inf_",i), boxplot_fct_inf_status)
}

boxplot_inf_variant <- plot_grid(
  boxplot_inf_Alpha,
  boxplot_inf_Delta,
  `boxplot_inf_Omicron BA.1`,
  `boxplot_inf_Omicron BA.2`,
  `boxplot_inf_Omicron BA.4/5`,
  ncol = 1
)

boxplot_inf_variant %>%
  write_rds(file = "output/boxplot_inf_variant.rds")


# wGLM survey -------------------------------------------------------------

tab_inf_variant  <- tab_inf_variant  %>%
  mutate(svyglm = map2(.x = data_ips,
                       .y = variant,
                      .f = ~fun_svyglm(data_ips = .x,
                                       variant = .y)))




fun_plot_est <- function(svymod, variant) {
  plot_outcomes_fct_previous_infection <-  svymod %>%
    map(.f = \(x) {
      form <- x %>% tidytable::pull(formula) %>% as.character() %>% pluck(2) %>%
        str_remove_all("`")
      x %>% tidy(conf.int = T) %>%
        filter(str_detect(term, "fct_previous_infection")) %>%
        add_row(term = "no prior infection") %>%
        mutate(
          term = str_remove(term, "fct_previous_infection") %>%
            # voeg enter toe
            str_replace(pattern = "prior ", replacement = "prior\n") %>% 
            factor(levels = rev(c("no prior\ninfection", 
                              "1 prior\ninfection", 
                              "2+ prior\ninfections"))),
          label = if_else(term == "no prior\ninfection", "reference", NA)
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
        geom_text(x = 3) +
        scale_x_continuous(breaks = seq(-1, 10, 1)) +
        coord_cartesian(xlim = c(-1, 10)) +
        labs(y = "") +
        theme(legend.position = "none", 
              axis.text.y = element_text(size = 10),
              axis.title.x = element_blank(),
              plot.margin = unit(c(0,0.1,0,0), "cm"),
              plot.title = element_text(face = "plain", size = 10,
                                        margin = margin(t = 10, b = 0))) +
        ggtitle(label = variant)
    })
}

# plot result -------------------------------------------------------------

tab_inf_variant <- tab_inf_variant %>%
  mutate(plot = map2(
    .x = svyglm,
    .y = variant,
    .f = ~ fun_plot_est(.x, .y)
  ) %>% unlist(F, F))



# add tidy outcome --------------------------------------------------------


fun_tidy_outcome <- function(svymod) {
  tidy_outcome <-  svymod %>%
    map(.f = \(x) {
      x %>% tidy(conf.int = T) 
    })
}

tab_inf_variant <- tab_inf_variant %>%
  mutate(tidy = map(
    .x = svyglm,
    .f = ~ fun_tidy_outcome(.x)
  ) %>% unlist(F, F))


plot_prev_inf_variant <- tab_inf_variant %>% 
  pull(plot) %>%
  cowplot::plot_grid(plotlist = ., nrow = 5)

# save --------------------------------------------------------------------
tab_inf_variant %>% saveRDS("cache/tab_inf_variant.rds")

plot_prev_inf_variant %>% 
  write_rds(file = "output/plot_prev_inf_variant.rds")

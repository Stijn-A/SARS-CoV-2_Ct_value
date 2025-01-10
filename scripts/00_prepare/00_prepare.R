lapply(
  c(
    "tidyverse",
    "tidymodels",
    "tidytable",
    "rlang",
    "tableone",
    "xlsx",
    "lubridate",
    "readxl",
    "cowplot",
    "vctrs",
    "mgcv",
    "survey",
    "splines",
    "tidytable"
  ),
  require,
  character.only = TRUE
)
theme_set(theme_bw())

dim_tab_treatment <- tibble(
  num_treatment =
    seq(0, 9, 1),
  fct_treatment =
    c(
      "unvacc::no prev inf",
      "unvacc::prev inf",
      "full vacc:<=90d:no prev inf",
      "full vacc:<=90d:prev inf",
      "full vacc:>90d:no prev inf",
      "full vacc:>90d:prev inf",
      "booster:<=90d:no prev inf",
      "booster:<=90d:prev inf",
      "booster:>90d:no prev inf",
      "booster:>90d:prev inf"
    ) %>% factor(., levels = .)
)

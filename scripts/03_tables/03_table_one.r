data_tab1 <- data_ct %>%
  select(
    lab,
    num_lab,
    `Ct ORF1ab`,
    `Ct N`,
    `Ct S`,
    age_group10,
    sex,
    testing_date,
    fct_vacc_status,
    fct_previous_infection,
    variant_short,
    variant_typing,
    test_month
  ) %>%
  mutate(variant_short = variant_short %>% fct_na_value_to_level("unknown"))
table1 <- bind_rows(
  data_tab1 %>%
    count(num_lab) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "total",
      variable = "",
      .before = 1
    ),
  # Ct values ---------------------------------------------------------------
  data_tab1 %>%
    group_by(num_lab) %>%
    summarise(
      min = min(`Ct N`, na.rm = TRUE),
      Q1 = quantile(`Ct N`, probs = 0.25, na.rm = TRUE),
      median = quantile(`Ct N`, probs = 0.5, na.rm = TRUE),
      Q3 = quantile(`Ct N`, probs = 0.75, na.rm = TRUE),
      max = max(`Ct N`, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(across(
      .col = !c(num_lab),
      .fns = \(.x) {
        round(.x, digits = 1)
      }
    )) %>%
    pivot_longer(cols = !c(num_lab),
                 names_to = "variable") %>%
    pivot_wider(names_from = c(num_lab),
                values_from = value) %>%
    mutate(
      across(.fns = as.character),
      type = "Ct N",
      variable = variable %>% factor(levels = c("min", "Q1", "median", "Q3",
                                                "max")),
      .before = 1
    ) %>%
    arrange(variable),
  # age group ---------------------------------------------------------------
  data_tab1 %>%
    count(num_lab, age_group10) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "Age group",
      .before = 1
    ) %>%
    rename(variable = age_group10) %>%
    arrange(variable),
  # sex ---------------------------------------------------------------------
  data_tab1 %>%
    count(num_lab, sex) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "sex",
      .before = 1
    ) %>%
    rename(variable = sex),
  # min max testing date ----------------------------------------------------
  data_tab1 %>%
    group_by(num_lab) %>%
    summarise(
      min_testing_date = min(testing_date),
      max_testing_date = max(testing_date)
    ) %>%
    pivot_longer(
      cols = c(min_testing_date, max_testing_date),
      names_to = "variable"
    ) %>%
    pivot_wider(names_from = c(num_lab),
                values_from = value) %>%
    mutate(
      across(.fns = as.character),
      type = "testing date",
      variable = variable %>% factor(levels = c("min_testing_date",
                                                "max_testing_date")),
      .before = 1
    ) %>%
    arrange(variable),
  # vaccination status ---------------------------------------------------------
  data_tab1 %>%
    count(num_lab, fct_vacc_status) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "fct_vacc_status",
      .before = 1
    ) %>%
    rename(variable = fct_vacc_status) %>%
    arrange(variable),
  # previous infections status ----------------------------------------------
  data_tab1 %>%
    count(num_lab, fct_previous_infection) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "fct_previous_infection",
      .before = 1
    ) %>%
    rename(variable = fct_previous_infection) %>%
    arrange(variable),
  # variant -----------------------------------------------------------------
  data_tab1 %>%
    count(num_lab, variant_short) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "variant_short",
      .before = 1
    ) %>%
    rename(variable = variant_short) %>%
    arrange(variable),
  # variant typing ----------------------------------------------------------
  data_tab1 %>%
    count(num_lab, variant_typing) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      across(.fns = as.character),
      type = "variant_typing",
      .before = 1
    ) %>%
    rename(variable = variant_typing) %>%
    arrange(variable),
) %>%
  group_by(type) %>%
  mutate(across(
    .cols = !c(variable),
    .fns = \(x) {
      round(as.numeric(x) / sum(as.numeric(x)) * 100, 1)
    },
    .names = "{.col}_perc"
  )) %>%
  ungroup
# order colnames
table1 <-
  table1 %>% relocate(type, variable, order(colnames(table1)))

table1 <- table1 %>%
  mutate(
    Lab_1_N_c = paste0(`Lab 1_N`, " (", `Lab 1_N_perc`, ")"),
    Lab_1_ORF1ab_c = paste0(`Lab 1_ORF1ab`, " (", `Lab 1_ORF1ab_perc`, ")"),
    Lab_1_S_c = paste0(`Lab 1_S`, " (", `Lab 1_S_perc`, ")"),
    Lab_2_N_c = paste0(`Lab 2_N`, " (", `Lab 2_N_perc`, ")"),
    Lab_2_ORF1ab_c = paste0(`Lab 2_ORF1ab`, " (", `Lab 2_ORF1ab_perc`, ")"),
    Lab_2_S_c = paste0(`Lab 2_S`, " (", `Lab 2_S_perc`, ")"),
    across(.cols = starts_with("Lab"),
           ~ str_remove(.x, pattern = "\\(NA\\)"))
  ) %>%
  select(type, variable, ends_with("_c"))
# supplement table --------------------------------------------------------
sup_table_one <- bind_rows(
  # test month --------------------------------------------------------------
  data_tab1 %>%
    count(num_lab, test_month) %>%
    pivot_wider(
      names_from = c(num_lab),
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(type = "test_month", .before = 1) %>%
    rename(variable = test_month),
)

#clean data
source("scripts/02_transform/koppeltabel_variant_tijd.R")

data_sgtf <- data_sgtf_org %>%
  filter(
    `S result` != "Negative",
    #remove missing Ct values
    !(is.na(`Ct N`) & is.na(`Ct ORF1ab`) & is.na(`Ct S`))
  ) %>%
  select(Monsternummer, `Ct ORF1ab`, `Ct S`, `Ct N`, `S result`, lab = Lab)

data_ct <- data_sgtf %>%
  # link epi data for each test
  left_join(
    data_teststr_org,
    by = "Monsternummer"
  ) %>%
  filter(Uitslag == "Positief")

logger::log_info(str_c(
  "Nr of positive tests ",
  data_ct %>%
    filter(Afspraak_start_datum %in% seq(as_date("2021-03-01"),
                                         as_date("2022-12-31"),
                                         by = "day")) %>%
    nrow
))

data_teststr <-  data_teststr_org %>%
  filter(Uitslag == "Positief" &
           !is.na(Pseudoniem) & Pseudoniem %in% data_ct$Pseudoniem) %>%
  select(Pseudoniem, Afspraak_start_datum, Monsternummer)

data_teststr <-  data_teststr %>%
  arrange(Pseudoniem, Afspraak_start_datum) %>%
  mutate(interval_afspraken =
            as.numeric(Afspraak_start_datum - lag(Afspraak_start_datum)),
          .by = Pseudoniem) %>%
  filter(interval_afspraken > 29 | is.na(interval_afspraken)) %>%
  mutate(
    n_pos_tests = row_number(),
    lag_monsternummer = lag(Monsternummer),
    lag_afspraak_start_datum = lag(Afspraak_start_datum),
    min_afspraak_start_datum = min(Afspraak_start_datum),
    .by = Pseudoniem
  )

data_teststr <- data_teststr %>%
  # remove 2 double monsternummers
  slice(1, .by = Monsternummer)

# data enrichment
data_ct <- data_ct %>%
  # remove monsternummer without any epi information
  filter(Monsternummer %in% data_teststr$Monsternummer) %>%
  left_join(
    data_teststr %>% select(!c(Pseudoniem, Afspraak_start_datum)),
    by = "Monsternummer"
  ) %>%
  left_join(
    koppeltabel_variant_S_result %>% select(
      testing_date = date,
      `S result`,
      variant_sgtf
    ),
    by = c("Afspraak_start_datum" = "testing_date", "S result")
  ) %>%
  left_join(
    data_kiemsurv %>%
      filter(!is.na(`CoronIT Id`)) %>%
      select(`CoronIT Id`,
             `VOC/VOI`,
             VOC,
             Clade),
    by = c("Monsternummer" = "CoronIT Id")
  )

# data enrichment variant of prior infection
data_ct <- data_ct %>%
  left_join(
    data_ct %>%
      select(Monsternummer,
             lag_variant_sgtf = variant_sgtf,
             lag_S_result = `S result`),
    by = c("lag_monsternummer" = "Monsternummer")
  ) %>%
  left_join(
    data_kiemsurv %>%
      filter(!is.na(`CoronIT Id`)) %>%
      select(`CoronIT Id`,
             lag_voc = VOC),
    by = c("lag_monsternummer" = "CoronIT Id")
  ) %>%
  mutate(
    lag_variant = coalesce(lag_variant_sgtf, lag_voc) %>%
      factor(
        levels = c(
          "Pre-VOC",
          "Alpha",
          "Beta",
          "Gamma",
          "Delta",
          "Omicron BA.1",
          "Omicron BA.2",
          "Omicron BA.4",
          "Omicron BA.5",
          "Omicron BA.4/5",
          "Omicron BA.2.12.1",
          "Omicron BA.2.75",
          "Omicron BQ.1.1",
          "Omicron XBB.1"
        )
      ),
    lag_variant_short = case_when(
      lag_variant %in% c("BA.2",
                         "Omicron BA.2.12.1",
                         "Omicron BA.2.75") ~ "BA.2",
      lag_variant %in% c("BA.4", "BA.5", "BA.4/5") ~ "BA.4/5",
      TRUE ~ lag_variant
    ) %>%
      factor(
        levels = c(
          "Pre-VOC",
          "Alpha",
          "Beta",
          "Gamma",
          "Delta",
          "Omicron BA.1",
          "Omicron BA.2",
          "Omicron BA.4/5",
          "Omicron BQ.1.1",
          "Omicron XBB.1"
        )
      )
  ) %>%
  select(-c(lag_voc, lag_variant))

data_ct <- data_ct %>%
  mutate(
    testing_date = Afspraak_start_datum,
    test_month = floor_date(testing_date, "month") %>% factor,
    test_week = floor_date(testing_date, "week", week_start = 1),
    num_testing_date = as.numeric(Afspraak_start_datum),
    num_week = case_when(
      isoyear(testing_date) == 2021 ~ isoweek(testing_date),
      isoyear(testing_date) == 2022 ~
        isoweek(testing_date) + isoweek("2021-12-31"),
      isoyear(testing_date) == 2023 ~
        isoweek(testing_date) + isoweek("2022-12-31")
    ),
    age_group10 = Age %>% cut(
      breaks = c(18, 30, 40, 50, 60, 70, 80, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c(
        "18-29",
        "30-39",
        "40-49",
        "50-59",
        "60-69",
        "70-79",
        "80+"
      )
    ) %>%
      factor(levels = c(
        "18-29",
        "30-39",
        "40-49",
        "50-59",
        "60-69",
        "70-79",
        "80+"
      )),
    age_group10_ref = age_group10 %>% factor(
      levels = c(
        "40-49",
        "18-29",
        "30-39",
        "50-59",
        "60-69",
        "70-79",
        "80+"
      )
    ),
    # if WGS variant is not availible sgtf variant
    variant = coalesce(VOC, variant_sgtf) %>%
      factor(levels = c("Pre-VOC", "Alpha", "Beta", "Gamma", "Delta",
                        "Omicron BA.1", "Omicron BA.2", "Omicron BA.4",
                        "Omicron BA.5", "Omicron BA.4/5",  "Omicron BA.2.12.1",
                        "Omicron BA.2.75", "Omicron BQ.1.1", "Omicron XBB.1")),
    variant_short = case_when(
      variant %in% c("Omicron BA.2", "Omicron BA.2.12.1", "Omicron BA.2.75") ~
        "Omicron BA.2",
      variant %in% c("Omicron BA.4", "Omicron BA.5", "Omicron BA.4/5") ~
        "Omicron BA.4/5",
      TRUE ~ variant
    ) %>%
      factor(levels = c("Pre-VOC", "Alpha", "Beta", "Gamma", "Delta",
                        "Omicron BA.1", "Omicron BA.2", "Omicron BA.4/5",
                        "Omicron BQ.1.1", "Omicron XBB.1")),
    variant_short_ref = variant_short %>%
      factor(levels = c("Delta", "Pre-VOC", "Alpha", "Beta", "Gamma",
                        "Omicron BA.1", "Omicron BA.2", "Omicron BA.4/5",
                        "Omicron BQ.1.1", "Omicron XBB.1")),
    variant_typing = case_when(
      !is.na(variant_short) & !is.na(VOC) ~
        "WGS", !is.na(variant_short) &
        !is.na(variant_sgtf) ~
        "SGTF",
      is.na(variant_short) ~
        "not typed"
    ),
    time_since_symp_onset = if_else(Onset_date == DUFS_EZD_0_42,
                                    as.numeric(testing_date - DUFS_EZD_0_42),
                                    NaN),
    time_since_symp_onset_cat = time_since_symp_onset %>% cut(
      breaks = c(0, 3, 6, 9, 13, 19, 25, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c("0-2", "3-5", "6-8", "9-12", "13-18", "19-24", "25-42")
    ) %>%
      factor(levels = c("0-2", "3-5", "6-8", "9-12",
                        "13-18", "19-24", "25-42")),
    time_since_previous_infection = interval_afspraken %>% cut(
      breaks = c(30, 90, 150, 240, 450, 750, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c("30-89", "90-149", "150-239", "240-449", "450-749", "750+")
    ) %>%
      fct_na_value_to_level("No prior\ninfection") %>%
      factor(levels = c("30-89", "90-149", "150-239", "240-449", "450-749",
                        "750+", "No prior\ninfection")),
    interval_vaccination =
      as.numeric(Afspraak_start_datum - Vaccination_date_latest),
    time_since_vaccination = interval_vaccination %>% cut(
      breaks = c(7, 60, 120, 180, 300, 600, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = c("7-59", "60-119", "120-179", "180-299", "300-599", "600+")
    ) %>%
    fct_na_value_to_level("Unvaccinated/Unknown") %>%
    factor(levels = c("7-59", "60-119", "120-179", "180-299", "300-599",
                      "600+", "Unvaccinated/Unknown")),
    Previous_infection = if_else(!is.na(lag_monsternummer),
                                 "Prior infection",
                                 "No prior infection") %>% as_factor(),
    n_previous_infections = if_else(
      n_pos_tests - 1 <= 2, n_pos_tests - 1, 3 # cat 0,1,2,3 (== 3+)
    ) %>% as.factor(),
    fct_previous_infection = case_when(
      n_pos_tests == 1 ~ "no prior infection",
      n_pos_tests == 2 ~ "1 prior infection",
      n_pos_tests > 2 ~ "2+ prior infections"
    ) %>% factor(levels = c("no prior infection",
                            "1 prior infection",
                            "2+ prior infections")),
    num_previous_infection = fct_previous_infection %>%
    case_match("no prior infection" ~ 0,
               "1 prior infection" ~ 1,
               "2+ prior infections" ~ 2),
    # use factor
    fct_vacc_status = case_match(
      Vaccination_status,
      "Ongevaccineerd" ~ "unvacc",
      "Volledig" ~ "full_vacc",
      "Booster" ~ "booster",
      TRUE ~ NA
    ) %>%
    factor(levels = c("unvacc", "full_vacc", "booster")),
    num_vacc_status = fct_vacc_status %>% case_match(
      "unvacc" ~ 0,
      "full_vacc" ~ 1,
      "booster" ~ 2
    ),
    fct_immune_status = interaction(fct_vacc_status, fct_previous_infection),
    fct_vacc_time = interaction(fct_vacc_status, time_since_vaccination) %>%
    recode(
      "unvacc.7-59" = NA,
      "unvacc.60-119" = NA,
      "unvacc.120-179" = NA,
      "unvacc.180-299" = NA,
      "unvacc.300-599" = NA,
      "full_vacc.Unvaccinated/Unknown" = NA,
      "unvacc.Unvaccinated/Unknown" = "unvacc",
    ) %>%
    str_replace(pattern = "\\.", replacement = " ") %>%
    factor(levels = c(
      "unvacc",
      "full_vacc 7-59",
      "full_vacc 60-119",
      "full_vacc 120-179",
      "full_vacc 180-299",
      "full_vacc 300-599",
      #level not present in data: "full_vacc.600+",
      "booster 7-59",
      "booster 60-119",
      "booster 120-179",
      "booster 180-299"#,
      #level not present in data: "booster.300-599"
    )),
    num_vacc_time = as.numeric(fct_vacc_time) - 1,
    prev_inf_time = interaction(fct_previous_infection,
                                time_since_previous_infection),
    # rename and combine groups
    fct_prev_inf_time = case_match(
      prev_inf_time,
      "no prior infection.No prior\ninfection" ~ "no prior infection",
      c("1 prior infection.450-749",
        "1 prior infection.750+") ~ "1 prior infection.450-849",
      c(
        "2+ prior infections.30-89",
        "2+ prior infections.90-149",
        "2+ prior infections.150-239"
      ) ~ "2+ prior infections.30-239",
      c(
        "2+ prior infections.240-449",
        "2+ prior infections.450-749"
      ) ~ "2+ prior infections.240-749",
      .default = prev_inf_time
    ) %>%
    str_replace(pattern = "\\.", replacement = "\n") %>%
    factor(
      levels = c(
        "no prior infection",
        "1 prior infection\n30-89",
        "1 prior infection\n90-149",
        "1 prior infection\n150-239",
        "1 prior infection\n240-449",
        "1 prior infection\n450-849",
        "2+ prior infections\n30-239",
        "2+ prior infections\n240-749"
      )
    ),
    num_prev_inf_time = as.numeric(fct_prev_inf_time) - 1,
    fct_prev_inf = Previous_infection %>% fct_recode(
      `prev inf` = "Prior infection",
      `no prev inf`  = "No prior infection"
    ),
    
    adjustment_prior_inf = case_when(
      fct_vacc_status == "unvacc" | fct_previous_infection == "no prior infection" ~ fct_previous_infection,
      lag_afspraak_start_datum < Vaccination_date_latest ~ fct_previous_infection, # if before vacc, adjust for confounding
      fct_previous_infection == "1 prior infection" ~ "no prior infection", # one inf after vacc, no prior inf
      min_afspraak_start_datum < Vaccination_date_latest ~ "1 prior infection",
      !min_afspraak_start_datum < Vaccination_date_latest ~ "no prior infection"
      ),
    
    adjustment_vacc = case_when(
      fct_vacc_status == "unvacc" | fct_previous_infection == "no prior infection" ~ fct_vacc_status,
      lag_afspraak_start_datum > Vaccination_date_latest ~ fct_vacc_status, # when vaccination was before prior infection
      !lag_afspraak_start_datum > Vaccination_date_latest & fct_vacc_status == "full_vacc" ~ "unvacc",
      !lag_afspraak_start_datum > Vaccination_date_latest & fct_vacc_status == "booster" &
        (lag_afspraak_start_datum - Vaccination_date_latest) < 365 ~ "full_vacc",
      TRUE ~ "unvacc"
    ),
    
    target = case_when(
      lab %in% c("Saltro", "Synlab") & !is.na(`Ct S`) ~ "S,N,ORF1ab",
      lab %in% c("Saltro", "Synlab") & is.na(`Ct S`) ~ "N,ORF1ab"
    ),
    lab = lab %>% factor(),
    num_lab = case_match(
      lab,
      "Saltro" ~ "Lab 1",
      "Synlab" ~ "Lab 2"
    ),
    sex = sex %>% case_match(
      "Man" ~ "male",
      "Vrouw" ~ "female",
      "Niet vermeld" ~ "unknown"
    )
  )

data_ct <- data_ct %>%
  filter(test_month %in% seq(as_date("2021-03-01"),
                             as_date("2022-12-01"),
                             by = "month"))

source("scripts/04_figures/04_flowchart_exclusions.R")

# selection
data_ct <- data_ct %>% filter(
  Age >= 18,
  sex != "unknown",
  Symptoms_at_least_1,
  time_since_symp_onset <= 24,
  !is.na(variant_short),
  # with known vaccination status and time since vaccination
  !is.na(num_vacc_time),
  round(`Ct ORF1ab`) <= 30 & round(`Ct N`) <= 30,
  # remove recently (<7 days) boostered from full vacc
  !(fct_vacc_status == "full_vacc" & Vaccination_nr_dose == 3)
) %>%
  rename(age = Age) %>% 
  droplevels()

data_ct <- data_ct %>%
  select(!Pseudoniem)
# To-do add variant prior infection

logger::log_info(str_c(
  "Nr of positive tests with epi and variant information ",
  data_ct %>%
    nrow
))

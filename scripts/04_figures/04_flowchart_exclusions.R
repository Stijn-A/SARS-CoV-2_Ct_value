library(flowchart)

# keep track of nr of filtered samples per step
tab_exlusion_list <-
  list(
    list(filter = "pos_test", value = data_ct %>% nrow()),
    list(
      filter = "excluded_age",
      value = data_ct %>% filter(Age >= 18) %>% nrow()
    ),
    list(
      filter = "sex_unknown",
      value = data_ct %>% filter(Age >= 18,
                                 sex != "unknown") %>% nrow()
    ),
    list(
      filter = "no symptoms",
      value = data_ct %>% filter(Age >= 18,
                                 sex != "unknown",
                                 Symptoms_at_least_1) %>% nrow()
    ),
    list(
      filter = "time_since_symptom_onset",
      value = data_ct %>% filter(
        Age >= 18,
        sex != "unknown",
        Symptoms_at_least_1,
        time_since_symp_onset <= 24
      ) %>% nrow()
    ),
    list(
      filter = "missing variant",
      value = data_ct %>% filter(
        Age >= 18,
        sex != "unknown",
        Symptoms_at_least_1,
        time_since_symp_onset <= 24,
        !is.na(variant_short)) %>% nrow()),
    list(
      filter = "missing vaccination data",
      value = data_ct %>% filter(
        Age >= 18,
        sex != "unknown",
        Symptoms_at_least_1,
        time_since_symp_onset <= 24,
        !is.na(variant_short),
        # with known vaccination status and time since vaccination
        !is.na(num_vacc_time)) %>% nrow()),
    list(filter = "Ct above 30",
         value = data_ct %>% filter(
           Age >= 18,
           sex != "unknown",
           Symptoms_at_least_1,
           time_since_symp_onset <= 24,
           !is.na(variant_short),
           # with known vaccination status and time since vaccination
           !is.na(num_vacc_time),
           round(`Ct ORF1ab`) <= 30 & round(`Ct N`) <= 30) %>% nrow()),
    list(
      filter = "booster_within_7days",
      value = data_ct %>% filter(
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
      ) %>% nrow()
    ))

tab_exlusion <- bind_rows(lapply(tab_exlusion_list, as.data.frame)) %>% 
  mutate(
    removed = value - lag(value)
  )


flowchart_exclusions <- data_ct %>% 
  as_fc(label = "Positive tests") %>% 
  fc_filter(Age >= 18, label = "Age >= 18 years", show_exc = TRUE) %>% 
  fc_filter(sex != "unknown", label = "Data on sex\nnot missing", show_exc = TRUE) %>% 
  fc_filter(Symptoms_at_least_1, label = "Symptomatic infection", show_exc = TRUE) %>% 
  fc_filter(time_since_symp_onset <= 24, label = "Sample with known symptom onset and\ntaken within 24 days after\nsymptom onset", show_exc = TRUE) %>% 
  fc_filter(!is.na(variant_short), label = "Data on variant\nnot missing", show_exc = TRUE) %>% 
  fc_filter(!is.na(num_vacc_time), label = "Data on vaccination\nnot missing", show_exc = TRUE) %>% 
  fc_filter(round(`Ct ORF1ab`) <= 30 & round(`Ct N`) <= 30, label = "Tests with a signal of Ct value\n≤ 30 for ORF1ab and N", show_exc = TRUE) %>% 
  fc_filter(!(fct_vacc_status == "full_vacc" & Vaccination_nr_dose == 3), label = "Without inconclusive vaccination status", show_exc = TRUE) %>% 
  fc_draw()

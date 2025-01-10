mod_tidy <- function(lm_tidy, str_formula) {
  
  variable = gsub(".*\\~ ", "", str_formula)
  
  lm_tidy = lm_tidy %>% 
    filter(term != "(Intercept)")
  
  if (variable == "age_group10_ref") {
    lm_tidy2 = lm_tidy %>%
      mutate(
        ylab_term = term %>% str_remove("age_group10_ref") %>%
          str_remove_all("`") %>%
          factor(
            levels = c(
              "18-29",
              "30-39","40-49",
              "50-59",
              "60-69",
              "70-79",
              "80+"
            )
          )
      )
  } else if (variable == "sex") {
    lm_tidy2 = lm_tidy %>%
      mutate(
        ylab_term = term %>% str_remove("sex") %>%
          
          factor(
            levels = c("female", "male")
          )
      )
  } else if (variable == "num_lab") {
    lm_tidy2 = lm_tidy %>%
      mutate(
        ylab_term = term %>% str_remove("num_lab") %>% str_remove_all("`") %>%
          factor(levels = data_ct_org$num_lab %>% unique())
      )
  } else if (variable == "time_since_symp_onset_cat") {
    lm_tidy2 = lm_tidy %>%
      mutate(
        ylab_term = term %>% str_remove("time_since_symp_onset_cat") %>%
          str_remove_all("`") %>%
          factor(
            levels = c("0-2", "3-5", "6-8", "9-12", "13-18", "19-24")
          )
      )
  } else if (variable == "fct_vacc_status") {
    lm_tidy2 = lm_tidy %>%
      mutate(
        ylab_term = term %>% str_remove("fct_vacc_status") %>%
          factor(levels = data_ct_org$fct_vacc_status %>% levels)
      )
  } else if (variable == "variant_short_ref") {
    lm_tidy2 = lm_tidy %>%
      mutate(ylab_term = term %>% str_remove("variant_short_ref") %>%
               str_remove_all("`") %>%
               factor(
                 levels = data_ct_org$variant_short %>% levels
               ))
  } else if (variable == "fct_prev_inf_time") {
    lm_tidy2 = lm_tidy %>%
      mutate(ylab_term = term %>% str_remove("fct_prev_inf_time") %>%
               str_remove_all("`") %>%
               gsub(pattern = "[\\\n]",replacement =  "test" ) %>% 
               str_replace("testn", "\n") %>% 
               factor(
                 levels = data_ct_org$fct_prev_inf_time %>% levels
               )
             )
  } else if (variable == "fct_vacc_time") {
    lm_tidy2 = lm_tidy %>%
      mutate(ylab_term = term %>% str_remove("fct_vacc_time") %>%
               str_remove_all("`") %>%
               factor(
                 levels = data_ct_org$fct_vacc_time %>% levels
               ))
    
  } else {
    lm_tidy2 = lm_tidy %>%
      mutate(ylab_term = term %>% str_remove(variable))
  }
  return(lm_tidy2)
}

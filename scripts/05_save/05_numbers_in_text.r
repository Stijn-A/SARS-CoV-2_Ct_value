#
# Results -----------------------------------------------------------------
#

# Study population --------------------------------------------------------
data_ct %>% nrow
data_ct %>% 
  count(`S result`) %>% 
  mutate(
    `%` = n / sum(n) * 100
  )

# Exploratory analysis ----------------------------------------------------


# Vaccination status  -----------------------------------------------------
c("cache/mod_outcome_vacc.rds", "cache/mod_outcome_vacc_time.rds",
  "cache/mod_outcome_inf.rds", "cache/mod_outcome_inf_time.rds",
  "cache/tab_inf_variant.rds", "cache/tab_vacc_variant.rds") %>% 
  map(.x = ., .f = ~tail(file.info(.x)$ctime))

mod_outcome_vacc <- readRDS("cache/mod_outcome_vacc.rds")
mod_outcome_vacc_time <- readRDS("cache/mod_outcome_vacc_time.rds")
mod_outcome_vacc %>% 
  filter(str_detect(term, "vacc_status"))

mod_outcome_vacc_time %>% 
  filter(str_detect(term, "vacc_time"))


# Previous infection status  ----------------------------------------------
mod_outcome_inf <- readRDS("cache/mod_outcome_inf.rds")
mod_outcome_inf_time <- readRDS("cache/mod_outcome_inf_time.rds")

mod_outcome_inf %>% 
  filter(str_detect(term, "previous_infection"))

mod_outcome_inf_time %>% 
  filter(str_detect(term, "prev_inf_time"))

# Effect stratified by variant of current infection -----------------------
tab_inf_variant <- readRDS("cache/tab_inf_variant.rds")

tab_inf_variant %>% 
  filter(variant == "Alpha") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))

tab_inf_variant %>% 
  filter(variant == "Delta") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.1") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.2") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.4/5") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


# vaccination & variant ---------------------------------------------------

tab_vacc_variant <- readRDS("cache/tab_vacc_variant.rds")

tab_vacc_variant %>% 
  filter(variant == "Alpha") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))

tab_inf_variant %>% 
  filter(variant == "Delta") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.1") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.2") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))


tab_inf_variant %>% 
  filter(variant == "Omicron BA.4/5") %>% 
  pull(tidy) %>% 
  select(term, estimate, conf.low, conf.high) %>% 
  filter(str_detect(term, pattern = "fct_previous_infection")) %>% 
  mutate(
    across(where(is.numeric), ~ round(., 2)))



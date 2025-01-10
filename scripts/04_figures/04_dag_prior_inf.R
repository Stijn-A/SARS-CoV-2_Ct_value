library(tidyverse)
library(dagitty)
library(ggdag)

# Specify a simple DAG containing one path
dag_prior_inf <- dagitty(
  "dag{
prior_infection -> Ct
vacc_before_prior_infection -> Ct
vacc_after_prior_infection -> Ct
age -> Ct
lab -> Ct
sex -> prior_infection
test_date -> Ct
time_onset -> Ct
variant -> Ct

prior_infection <- vacc_before_prior_infection
prior_infection -> vacc_after_prior_infection
vacc_before_prior_infection -> vacc_after_prior_infection

lab -> prior_infection
age -> prior_infection
test_date -> prior_infection
time_onset <- prior_infection
prior_infection -> variant
test_date -> prior_infection
variant <- vacc_before_prior_infection
variant <- vacc_after_prior_infection
                          }"
) %>%
  tidy_dagitty()

dag_prior_inf <- dag_prior_inf %>%
  as_tibble() %>%
  mutate(
    variable_type = case_when(
      name %in% c(
        "vacc_before_prior_infection",
        "age",
        "lab",
        "sex",
        "test_date",
        "time_onset"
      ) ~ "confounding variables",
      name %in% c(
        "vacc_after_prior_infection",
        "variant"
      ) ~ "mediator",
      
      name == "Ct" ~ "outcome",
      name == "prior_infection" ~ "exposure"
    ),
    # Not yet in use
    studied_effect = if_else(name == "vacc_status" &
                               to == "Ct", # exposure to outcome edge
                             TRUE,
                             FALSE)
  )

dag_prior_inf %>%
  ggplot(aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend
  )) +
  geom_dag_point(aes(col = variable_type)) +
  geom_dag_edges() +
  geom_dag_text(col = "black") +
  scale_color_brewer(palette = "Set2") +
  theme_dag() +
  scale_color_hue(breaks = c("parent", "child")) #  ignores NA in legend

# small dag with outcome, exposure and prior infection
dag_prior_inf <- dagitty(
  "dag{
prior_infection -> Ct
variant -> Ct
prior_infection -> variant
                          }"
) %>%
  tidy_dagitty() %>% as_tibble() %>%
  mutate(
    variable_type = case_when(
      name %in% c(
        "vacc_before_prior_infection",
        "age",
        "lab",
        "sex",
        "test_date",
        "time_onset"
      ) ~ "confounding variables",
      name %in% c(
        "vacc_after_prior_infection",
        "variant"
      ) ~ "mediator",
      
      name == "Ct" ~ "outcome",
      name == "prior_infection" ~ "exposure"
    ),
    # Not yet in use
    studied_effect = if_else(name == "vacc_status" &
                               to == "Ct", # exposure to outcome edge
                             TRUE,
                             FALSE)
  )

dag_prior_inf %>%
  ggplot(aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend
  )) +
  geom_dag_point(aes(col = variable_type)) +
  geom_dag_edges() +
  geom_dag_text(col = "black") +
  scale_color_brewer(palette = "Set2") +
  theme_dag() +
  scale_color_hue(breaks = c("parent", "child")) #  ignores NA in legend

library(tidyverse)
library(dagitty)
library(ggdag)

# Specify a simple DAG containing one path
dag_vacc <- dagitty(
  "dag{
vacc_status -> Ct
prior_status_before_vacc -> Ct
prior_status_after_vacc -> Ct
age -> Ct
lab -> Ct
sex -> vacc_status
test_date -> Ct
time_onset -> Ct
variant -> Ct

vacc_status <- prior_status_before_vacc
vacc_status -> prior_status_after_vacc
prior_status_before_vacc -> prior_status_after_vacc

lab -> vacc_status
age -> vacc_status
test_date -> vacc_status
time_onset <- vacc_status
vacc_status -> variant
test_date -> vacc_status
variant <- prior_status_before_vacc
variant <- prior_status_after_vacc
                          }"
) %>%
  tidy_dagitty()

dag_vacc <- dag_vacc %>%
  as_tibble() %>%
  mutate(
    variable_type = case_when(
      name %in% c(
        "prior_status_before_vacc",
        "age",
        "lab",
        "sex",
        "test_date",
        "time_onset"
      ) ~ "confounding variables",
      name %in% c(
        "prior_status_after_vacc",
        "variant"
      ) ~ "mediator",
      
      name == "Ct" ~ "outcome",
      name == "vacc_status" ~ "exposure"
    ),
    # Not yet in use
    studied_effect = if_else(name == "vacc_status" &
                               to == "Ct", # exposure to outcome edge
                             TRUE,
                             FALSE)
  )

dag_vacc %>%
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
dag_vacc <- dagitty(
  "dag{
vacc_status -> Ct
prior_inf_before_vacc -> Ct
prior_inf_after_vacc -> Ct

vacc_status <- prior_inf_before_vacc
vacc_status -> prior_inf_after_vacc
prior_status_before_vacc -> prior_inf_after_vacc
                          }"
) %>%
  tidy_dagitty() %>% 
  as_tibble() %>%
  mutate(
    variable_type = case_when(
      name %in% c(
        "prior_inf_before_vacc",
        "age",
        "lab",
        "sex",
        "test_date",
        "time_onset"
      ) ~ "confounding variables",
      name %in% c(
        "prior_status_after_vacc",
        "variant"
      ) ~ "mediator",
      
      name == "Ct" ~ "outcome",
      name == "vacc_status" ~ "exposure"
    ),
    # Not yet in use
    studied_effect = if_else(name == "vacc_status" &
                               to == "Ct", # exposure to outcome edge
                             TRUE,
                             FALSE)
  )

dag_vacc %>%
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

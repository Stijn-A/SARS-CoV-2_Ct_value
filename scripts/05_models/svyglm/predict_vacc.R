library(survey)
library(splines)

mod_outcome_vacc <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/cache/mod_outcome_vacc.rds")


vacc.preddata <- expand_grid(
  fct_vacc_status = factor(x = c("full_vacc", "booster"), 
                           levels = data_ct$fct_vacc_status %>% levels),
  age = 50,
  sex = "male",
  fct_previous_infection = factor("no prior infection",
                                  levels = data_ct$fct_previous_infection %>% levels),
  time_since_symp_onset = data_ct %>% pull(time_since_symp_onset) %>% mean(),
  lab = factor("Saltro",
               levels = data_ct$lab %>% levels),
  variant_short = factor(x = c("Delta", "Omicron BA.1"),
                         levels = data_ct$variant_short %>% levels),
  num_testing_date = 18996
)


tmp <- predict(
  object = models_targets[[1]],
  newdata = vacc.preddata,
  #re.form=~0,
  #type = "link",
  interval="confidence",
  limit=0.95) %>%
  as_tibble




library(splines)
#Fit part
fit.data <- data.frame(y=rnorm(30),x=rnorm(30))
fit.ns <- lm(y ~ ns(x,3),data=fit.data)

#Predict
pred.data <- data.frame(y=rnorm(10),x=rnorm(10))
pred.fit <- predict(fit.ns,interval="confidence",limit=0.95,data.frame(x=pred.data$x))

lr_lm_wf <- function(str_formula, data = ct_train) {
  lm_mod <- linear_reg() 
  
  lm_workflow <-
    workflow() %>%
    step_dummy(all_nominal_predictors()) %>%
    add_model(lm_mod) %>%
    add_formula(as.formula(str_formula)) %>%
    fit(data = data)
  
  logger::log_info(str_c("wf ", str_formula))
  return(lm_workflow)
}
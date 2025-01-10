tidy_ci <- function(lm_workflow, str_formula){
  lm_tidy <- lm_workflow %>%
    tidy() %>% 
    mutate(conf.low = estimate - (1.96 * std.error), conf.high = estimate + (1.96 * std.error))
  
  logger::log_info(str_c("tidy ", str_formula))
  return(lm_tidy)}

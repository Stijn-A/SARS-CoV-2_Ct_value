plot_estimates_multivar <- function(glm_tidy, str_formula, y_filter) {
  if (exists("y_filter")) {
    glm_tidy <- glm_tidy %>%
      filter(str_detect(ylab_term, y_filter)) %>%
      mutate(
        ylab_term = ylab_term %>% str_remove(y_filter) %>%
          str_remove_all("`")
        %>% factor(levels = data_ct$fct_prev_inf_time %>% levels)
      )
  }
  
  object <- ggplot(glm_tidy,
                   aes(
                     estimate,
                     ylab_term,
                     xmin = conf.low,
                     xmax = conf.high,
                     height = 0
                   )) +
    geom_point() +
    geom_vline(xintercept = 0, lty = 4) +
    geom_errorbarh() +
    scale_color_manual(values = c("TRUE" =  "green", "FALSE" = "red")) +
    scale_fill_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    ylab(NULL) +
    ggtitle(str_formula) +
    theme(legend.position = "none")
  
  logger::log_info(str_c("plot ", str_formula))
  return(object)
}

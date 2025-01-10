plot_estimates_univar <- function(lm_tidy, str_formula) {
  
  variable = gsub(".*\\~ ", "", str_formula)
  
  ct = gsub(" \\~ .*", "", str_formula)
  
  object <- ggplot(lm_tidy,
                   aes(
                     estimate,
                     ylab_term,
                     xmin = conf.low,
                     xmax = conf.high,
                     height = 0#,
                     #color = sign_bonferoni
                   )) +
    geom_point() +
    geom_vline(xintercept = 0, lty = 4) +
    geom_errorbarh() +
    scale_color_manual(values = c("TRUE" =  "green", "FALSE" = "red")) +
    scale_fill_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE) +
    scale_x_continuous(breaks = seq(-5,5,1)) +
    ylab(variable) +
    ggtitle(ct) +
    coord_cartesian(xlim = c(-5, 5)) +
    theme(legend.position="none")
  
  logger::log_info(str_c("plot ", str_formula))
  return(object)
}

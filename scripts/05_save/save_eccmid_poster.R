

ggsave(
  str_c(
    path_output,
    'figure_vacc_inf_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".pdf"
  ),
  plot = plot_vacc_inf,
 device = "pdf",
 units = "in",
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  str_c(
    path_output,
    'figure_covariates_N_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".pdf"
  ),
  plot = figure_covariates_N,
  device = "pdf",
  units = "in",
  width = 18,
  height = 12,
  dpi = 300
)


ggsave(
  str_c(
    path_output,
    'fig_variant_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".pdf"
  ),
  plot = plot_variants,
  device = "pdf",
  units = "in",
  width = 8,
  height = 8,
  dpi = 300
)

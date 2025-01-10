dir.create(
  "/rivm/n/andewegs/Documents/SARS-CoV-2/09_Ct_values_extreme_value_theory/output/",
  showWarnings = F
)

path_output <- "output/"
dir.create("output", showWarnings = F)
path_output <-
  "/rivm/n/andewegs/Documents/SARS-CoV-2/09_Ct_values_extreme_value_theory/output/"

plot_outcomes_fct_previous_infection <-
  readRDS(
    "output/plot_outcomes_fct_previous_infection_N.rds"
  )
plot_outcomes_fct_vacc_status <-
  readRDS(
    "output/plot_outcomes_fct_vacc_status_N.rds"
  )
plot_outcomes_fct_prev_inf_time <-
  readRDS(
    "output/plot_outcomes_fct_prev_inf_time_N.rds"
  )
# add vacc time
plot_outcomes_fct_vacc_time <-
  readRDS(
    "output/plot_outcomes_fct_vacc_time_N.rds"
  )

plot_vacc <-
  plot_grid(
    plot_outcomes_fct_vacc_status,
    NULL,
    plot_outcomes_fct_vacc_time,
    ncol = 1,
    rel_heights = c(0.3,-0.05, 0.7),
    #axis = "rl",
    align = "v",
    labels = c("A", " ", "C")
  )

plot_inf <-
  plot_grid(
    plot_outcomes_fct_previous_infection,
    NULL,
    plot_outcomes_fct_prev_inf_time,
    ncol = 1,
    rel_heights = c(0.3,-0.05, 0.7),
    align = "v",
    labels = c("B", " ", "D")
  )

plot_vacc_inf <-
  plot_grid(plot_vacc,
            plot_inf,
            align = "h",
            rel_widths = c(0.45, 0.56))

plot_vacc_inf <-
  ggdraw(
    add_sub(
      plot_vacc_inf,
      "Estimate",
      vpadding = grid::unit(0, "lines"),
      y = 5,
      x = 0.5,
      vjust = 4.5,
      size = 18
    )
  )

tiff(
  str_c(
    path_output,
    'figure_vacc_inf_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 10,
  height = 8,
  res = 300
)

plot_vacc_inf
dev.off()


tiff(
  str_c(
    path_output,
    'figure_covariates_N_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 18,
  height = 12,
  res = 300
)

figure_covariates_N
dev.off()


tiff(
  str_c(
    path_output,
    'figure_covariates_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 14,
  height = 14,
  res = 300
)

figure_covariates
dev.off()


#
# tiff(
#   str_c(
#     path_output,
#     'figure_IV_N_',
#     format(now(), format = "%Y%m%d_%H%M"),
#     ".tiff"
#   ),
#   units = "in",
#   width = 18,
#   height = 9,
#   res = 300
# )
#
# figure_IV_N
# dev.off()

# vaccination -------------------------------------------------------------


# previous infection ------------------------------------------------------



# vaccination variant -----------------------------------------------------
plot_prev_inf_variant <-
  readRDS("output/plot_prev_inf_variant.rds")

plot_vacc_variant <-
  readRDS("output/plot_vacc_variant.rds")

plot_variants <-
  plot_grid(
    plot_vacc_variant,
    # NULL,
    plot_prev_inf_variant,
    ncol = 2,
    rel_heights = c(0.5, #-0.05,
                    0.5),
    labels = c("A", "B"),
    align = "h",
    scale = 0.93
  ) +
  draw_label(
    "Estimate",
    x = 0.5,
    y =  0,
    vjust = -0.5,
    angle = 0
  )

tiff(
  str_c(
    path_output,
    'fig_variant_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 8,
  height = 8,
  res = 300
)

plot_variants %>% plot()
dev.off()


write.xlsx(
  table1,
  file = str_c(
    path_output,
    "tabel1_Ct_values_evt_characteristics_",
    format(now(), format = "%Y%m%d_%H%M"),
    ".xlsx"
  )
)
list(
  "vacc" =  tab_sensitivity_analysis_vacc,
  "vacc_time" = tab_sensitivity_analysis_vacc_time,
  "prev_inf" = tab_sensitivity_analysis_prev_inf,
  "prev_inf_time" = tab_sensitivity_analysis_prev_inf_time
) %>%
  openxlsx::write.xlsx(.,
                       file = str_c(
                         path_output,
                         "tab_sensitivity_analysis_",
                         format(now(), format = "%Y%m%d_%H%M"),
                         ".xlsx"
                       ))

#data_sgtf %>% saveRDS("/rivm/r/COVID-19/Toegang_externen/EPI_Ct_values_teststraten/data/data_sgtf_20230222.rds")



# covariates figures ------------------------------------------------------

univar_covariates <- readRDS("~/01_projects/02_sars_cov_2/ct_values_evt/output/plot_univar_covariates.rds")

tiff(
  str_c(
    path_output,
    'figure_univar_covariates_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 10,
  height = 16,
  res = 300
)

univar_covariates
dev.off()

# tiff(
#   str_c(
#     path_output,
#     'figure_univar_covariates_2_',
#     format(now(), format = "%Y%m%d_%H%M"),
#     ".tiff"
#   ),
#   units = "in",
#   width = 10,
#   height = 10,
#   res = 300
# )
# 
# univar_covariates_2
# dev.off()


tiff(
  str_c(
    path_output,
    'figure_plot_prev_inf_time_variant_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 8,
  height = 12,
  res = 300
)

plot_prev_inf_time_variant
dev.off()


write.xlsx(table_sensitivity,
           file = str_c(
             path_output,
             "table_sensitivity_",
             format(now(), format = "%Y%m%d_%H%M"),
             ".xlsx"
           ))

tiff(
  str_c(
    path_output,
    'figure_supplement_labs_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 10,
  height = 6,
  res = 300
)

bar_graph_labs
dev.off()

# result ORF1ab -----------------------------------------------------------

plot_outcomes_fct_previous_infection_ORF1ab <-
  readRDS(
    "output/plot_outcomes_fct_previous_infection_orf1ab.rds"
  )
plot_outcomes_fct_vacc_status_ORF1ab <-
  readRDS(
    "output/plot_outcomes_fct_vacc_status_orf1ab.rds"
  )
plot_outcomes_fct_prev_inf_time_ORF1ab <-
  readRDS(
    "output/plot_outcomes_fct_prev_inf_time_orf1ab.rds"
  )
# add vacc time
plot_outcomes_fct_vacc_time_ORF1ab <-
  readRDS(
    "output/plot_outcomes_fct_vacc_time_orf1ab.rds"
  )

plot_vacc_ORF1ab <-
  plot_grid(
    plot_outcomes_fct_vacc_status_ORF1ab,
    NULL,
    plot_outcomes_fct_vacc_time_ORF1ab,
    ncol = 1,
    rel_heights = c(0.3,-0.05, 0.7),
    #axis = "rl",
    align = "v",
    labels = c("A", " ", "C")
  )

plot_inf_ORF1ab <-
  plot_grid(
    plot_outcomes_fct_previous_infection_ORF1ab,
    NULL,
    plot_outcomes_fct_prev_inf_time_ORF1ab,
    ncol = 1,
    rel_heights = c(0.3,-0.05, 0.7),
    align = "v",
    labels = c("B", " ", "D")
  )

plot_vacc_inf_ORF1ab <-
  plot_grid(plot_vacc_ORF1ab,
            plot_inf_ORF1ab,
            align = "h",
            rel_widths = c(0.45, 0.56))

tiff(
  str_c(
    path_output,
    'figure_vacc_inf_ORF1ab_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 10,
  height = 8,
  res = 300
)

plot_vacc_inf_ORF1ab
dev.off()


# IPS ---------------------------------------------------------------------

boxplot_fct_vacc_status <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_fct_vacc_status.rds")
boxplot_vacc_time <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_vacc_time.rds")
boxplot_fct_previous_infection <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_fct_previous_infection.rds")
boxplot_prev_inf_time <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_prev_inf_time.rds")

fig_IPS <- plot_grid(
  boxplot_fct_vacc_status,
  boxplot_vacc_time,
  boxplot_fct_previous_infection,
  boxplot_prev_inf_time,
  labels = c("A", "B", "C", "D")
)

tiff(
  str_c(
    path_output,
    'figure_IPS_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 12,
  height = 12,
  res = 300
)

fig_IPS
dev.off()

boxplot_vacc_variant <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_vacc_variant.rds")
boxplot_inf_variant <- readRDS("~/01_projects/01_sars-cov-2/ct_values_evt/output/boxplot_inf_variant.rds")

fig_PS_variants <- plot_grid(
  boxplot_vacc_variant,
  boxplot_inf_variant,
  labels = c("A", "B")
)

tiff(
  str_c(
    path_output,
    'figure_PS_variants_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 12,
  height = 12,
  res = 300
)

fig_PS_variants
dev.off()

tiff(
  str_c(
    path_output,
    'figure_flowchart_exclusions_',
    format(now(), format = "%Y%m%d_%H%M"),
    ".tiff"
  ),
  units = "in",
  width = 12,
  height = 12,
  res = 300
)

flowchart_exclusions %>% fc_draw()
dev.off()

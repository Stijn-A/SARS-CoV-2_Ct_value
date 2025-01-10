# pivot longer ------------------------------------------------------------
data_sgtf_long <- data_ct %>%
  select(!c(`Ct S`, `Ct ORF1ab`)) %>%
  pivot_longer(
    cols = starts_with("Ct"),
    values_to = "Ct",
    names_to = "Target",
    names_prefix = "Ct "
  ) %>%
  mutate(
    Target = factor(Target),
    variant_short = variant_short %>% fct_expand("unknown"),
    variant_short = replace_na(variant_short, "unknown"),
    time_of_testing = case_when(
      test_month %in% seq(as_date("2021-03-01"),
                          as_date("2021-08-01 "),
                          by = "1 month") ~ "2021/03 -\n2021/08  ",
      test_month %in% seq(as_date("2021-09-01"),
                          as_date("2022-02-01 "),
                          by = "1 month") ~ "2021/09 -\n2022/02  ",
      test_month %in% seq(as_date("2022-03-01"),
                          as_date("2022-08-01 "),
                          by = "1 month") ~ "2022/03 -\n2022/08  ",
      test_month %in% seq(as_date("2022-09-01"),
                          as_date("2022-12-01 "),
                          by = "1 month") ~ "2022/09 -\n2022/12  ",
    ) %>% factor(
      levels = c(
        "2021/03 -\n2021/08  ",
        "2021/09 -\n2022/02  ",
        "2022/03 -\n2022/08  ",
        "2022/09 -\n2022/12  "
      )
    )
  )

# plot age (groups) -------------------------------------------------------

histogram_agegroup <- ggplot(data = data_sgtf_long,
                             aes(x = Ct, color = age_group10,
                                 fill = age_group10)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Blues",
                     direction = -1,
                     guide = "none") +
  scale_fill_brewer(palette = "Blues",
                    direction = -1,
                    guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_agegroup <- ggplot(data = data_sgtf_long,
                           aes(x = age_group10, y = Ct, fill = age_group10)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Blues",
                    guide = guide_legend(nrow = 2),
                    direction = -1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "age group") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_agegroup <-
  cowplot::plot_grid(
    histogram_agegroup,
    boxplot_agegroup,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_agegroup
rm(histogram_agegroup, boxplot_agegroup)



# filter on age >= 18 -----------------------------------------------------

data_sgtf_long <- data_sgtf_long %>%
  filter(age >= 18)

# plot sex ----------------------------------------------------------------


histogram_sex <- ggplot(data = data_sgtf_long,
                        aes(x = Ct, color = sex, fill = sex)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    position = "identity",
    alpha = 0.1
  ) +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A"), guide = "none") +
  scale_color_manual(values = c("#33A02C", "#B2DF8A"), guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(0, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_sex <- ggplot(data = data_sgtf_long,
                      aes(x = sex, y = Ct, fill = sex)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_manual(values = c("#33A02C", "#B2DF8A")) +
  facet_wrap(vars(Target), nrow = 1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_sex <-
  cowplot::plot_grid(
    histogram_sex,
    boxplot_sex,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_sex
rm(histogram_sex, boxplot_sex)

# plot lab ----------------------------------------------------------------


histogram_lab <- ggplot(data = data_sgtf_long,
                        aes(x = Ct, color = num_lab, fill = num_lab)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    position = "identity",
    alpha = 0.1
  ) +
  scale_fill_manual(values = c("#FB9A99", "#E31A1C"), guide = "none") +
  scale_color_manual(values = c("#FB9A99", "#E31A1C"), guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_lab <- ggplot(data = data_sgtf_long,
                      aes(x = num_lab, y = Ct, fill = num_lab)) +
  #geom_violin() +
  geom_boxplot(alpha = 0.9) +
  scale_fill_manual(values = c("#FB9A99", "#E31A1C")) +
  facet_wrap(vars(Target), nrow = 1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  xlab("laboratory") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_lab <-
  cowplot::plot_grid(
    histogram_lab,
    boxplot_lab,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_lab
rm(histogram_lab, boxplot_lab)

# test_time ---------------------------------------------------------------

histogram_time <- ggplot(data = data_sgtf_long,
                         aes(x = Ct, color = time_of_testing,
                             fill = time_of_testing)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Purples",
                     direction = -1,
                     guide = "none") +
  scale_fill_brewer(palette = "Purples",
                    direction = -1,
                    guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_time <- ggplot(data = data_sgtf_long,
                       aes(x = time_of_testing, y = Ct,
                           fill = time_of_testing)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Purples",
                    guide = guide_legend(nrow = 1),
                    direction = -1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "calendar time of testing") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_time <-
  cowplot::plot_grid(
    histogram_time,
    boxplot_time,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_time
rm(histogram_time, boxplot_time)


# variant -----------------------------------------------------------------

histogram_variant <- ggplot(data = data_sgtf_long,
                            aes(x = Ct, color = variant_short,
                                fill = variant_short)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Set3",
                     direction = -1,
                     guide = "none") +
  scale_fill_brewer(palette = "Set3",
                    direction = -1,
                    guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_variant <- ggplot(data = data_sgtf_long,
                          aes(x = variant_short, y = Ct,
                              fill = variant_short)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Set3",
                    guide = guide_legend(nrow = 2),
                    direction = -1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "variant") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_variant <-
  cowplot::plot_grid(
    histogram_variant,
    boxplot_variant,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_variant
rm(histogram_variant, boxplot_variant)


# time since onset --------------------------------------------------------

histogram_onset <- ggplot(data = data_sgtf_long,
                          aes(x = Ct, color = time_since_symp_onset_cat,
                              fill = time_since_symp_onset_cat)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "YlOrRd",
                     guide = "none",
                     direction = -1) +
  scale_fill_brewer(palette = "YlOrRd",
                    guide = "none",
                    direction = -1) +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_onset <- ggplot(data = data_sgtf_long,
                        aes(x = time_since_symp_onset_cat, y = Ct,
                            fill = time_since_symp_onset_cat)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "YlOrRd",
                    guide = guide_legend(nrow = 1),
                    direction = -1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "time since onset (days)") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )


plot_onset <-
  cowplot::plot_grid(
    histogram_onset,
    boxplot_onset,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_onset
rm(histogram_onset, boxplot_onset)


# previous infection ------------------------------------------------------
histogram_previous_infection <- ggplot(data = data_sgtf_long,
                                       aes(x = Ct,
                                           color = fct_previous_infection,
                                           fill = fct_previous_infection)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.9,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Greens",
                     guide = "none",
                     direction = 1) +
  scale_fill_brewer(palette = "Greens",
                    guide = "none",
                    direction = 1) +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_previous_infection <- ggplot(data = data_sgtf_long,
                                     aes(x = fct_previous_infection, y = Ct,
                                         fill = fct_previous_infection)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Greens",
                    guide = guide_legend(nrow = 1),
                    direction = 1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "previous infection status") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_previous_infection <-
  cowplot::plot_grid(
    histogram_previous_infection,
    boxplot_previous_infection,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_previous_infection
rm(histogram_previous_infection, boxplot_previous_infection)


# time since previous infection -------------------------------------------
histogram_previous_infection_time <- # nolint
  ggplot(
    data = data_sgtf_long %>% filter(!is.na(fct_prev_inf_time)),
    aes(x = Ct, color = fct_prev_inf_time, fill = fct_prev_inf_time)
  ) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Spectral",
                     guide = "none",
                     direction = -1) +
  scale_fill_brewer(palette = "Spectral",
                    guide = "none",
                    direction = -1) +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(1, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_previous_infection_time <- # nolint
  ggplot(data = data_sgtf_long %>% filter(!is.na(fct_prev_inf_time)),
         aes(x = fct_prev_inf_time, y = Ct, fill = fct_prev_inf_time)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Spectral",
                    guide = guide_legend(nrow = 2),
                    direction = -1) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "time since previous infection\n(days)") +
  coord_flip(ylim = c(1, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_previous_infection_time <-
  cowplot::plot_grid(
    histogram_previous_infection_time,
    boxplot_previous_infection_time,
    align = "v",
    nrow = 2,
    rel_heights = c(0.3, 0.7)
  )
plot_previous_infection_time
rm(histogram_previous_infection_time,
   boxplot_previous_infection_time)



# vaccination -------------------------------------------------------------

histogram_vaccination <- ggplot(data = data_sgtf_long,
                                aes(x = Ct, color = fct_vacc_status,
                                    fill = fct_vacc_status)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "Oranges", guide = "none") +
  scale_fill_brewer(palette = "Oranges", guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(0, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_vaccination <- ggplot(data = data_sgtf_long,
                              aes(x = fct_vacc_status, y = Ct,
                                  fill = fct_vacc_status)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "Oranges", guide = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "vaccination status") +
  coord_flip(ylim = c(0, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_vaccination <-
  cowplot::plot_grid(
    histogram_vaccination,
    boxplot_vaccination,
    align = "v",
    nrow = 2,
    rel_heights = c(0.4, 0.6)
  )
plot_vaccination
rm(histogram_vaccination, boxplot_vaccination)

# vaccination time -------------------------------------------------------------

histogram_vaccination_time <-
  ggplot(data = data_sgtf_long %>% filter(!is.na(fct_vacc_time)),
         aes(x = Ct, color = fct_vacc_time, fill = fct_vacc_time)) +
  geom_area(
    aes(y = after_stat(count)),
    binwidth = 1,
    stat = "bin",
    alpha = 0.1,
    position = "identity"
  ) +
  scale_color_brewer(palette = "RdYlGn", guide = "none") +
  scale_fill_brewer(palette = "RdYlGn", guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 7),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  coord_cartesian(xlim = c(0, 30)) +
  facet_wrap(vars(Target), nrow = 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 15)
  )

boxplot_vaccination_time <-
  ggplot(data = data_sgtf_long %>% filter(!is.na(fct_vacc_time)),
         aes(x = fct_vacc_time, y = Ct, fill = fct_vacc_time)) +
  geom_boxplot(alpha = 0.9)  +
  scale_fill_brewer(palette = "RdYlGn", guide = guide_legend(nrow = 2)) +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    minor_breaks = NULL,
    # decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05))
  ) +
  # decrease gap between x-axis and plot area (but keep padding on top)
  scale_x_discrete(expand = expansion(mult = c(0, .05))) +
  facet_wrap(vars(Target), nrow = 1) +
  labs(x = "time since vaccination\n(days)") +
  coord_flip(ylim = c(0, 30)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # remove facet labels
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, -10, -5, -10),
    text = element_text(size = 15)
  )

plot_vaccination_time <-
  cowplot::plot_grid(
    histogram_vaccination_time,
    boxplot_vaccination_time,
    align = "v",
    nrow = 2,
    rel_heights = c(0.3, 0.7)
  )
plot_vaccination_time
rm(histogram_vaccination_time, boxplot_vaccination_time)


# combine plots F1 -----------------------------------------------------------

figure_covariates_N <- cowplot::plot_grid(
  plot_agegroup,
  plot_sex,
  plot_lab,
  plot_onset,
  plot_variant,
  plot_time,
  labels = c("A", "B", "C",
             "D", "E", "F")
)


# combine plots F2 -----------------------------------------------------------


figure_IV_N <- cowplot::plot_grid(
  plot_vaccination,
  plot_previous_infection,
  plot_vaccination_time,
  plot_previous_infection_time,
  rel_heights = c(0.4, 0.6),
  labels = c("A", "B", "C", "D")
)

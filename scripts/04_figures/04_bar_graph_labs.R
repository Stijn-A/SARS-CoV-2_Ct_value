bar_graph_labs <- data_ct %>%
  # data wrangling
  mutate(Afspraak_start_week = floor_date(Afspraak_start_datum, unit = "week")) %>%
  count(Afspraak_start_week, num_lab) %>%
  complete(Afspraak_start_week, num_lab, fill = list(n = 0)) %>%
  # create basic ggplot object
  ggplot(aes(x = Afspraak_start_week,
             y = n,
             fill = num_lab)) +
  geom_col(# display bars side-to-side
    position = "dodge",) +
  # title and axes titles
  labs(x = "Testing date (weeks)",
       y = "Number of positive tests") +
  scale_x_date(date_breaks  = "2 months") +
  scale_fill_brewer(name = "",
                    type = "qual",
                    palette = "Set1") +
  # customize y-axis
  scale_y_continuous(# decrease gap between x-axis and plot area (but keep padding on top)
    expand = expansion(mult = c(0, .05)))

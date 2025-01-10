relational_table_clade_VOC <- tibble( # nolint
  Clade = c(
    "19A",
    "19B",
    "20A",
    "20B",
    "20C",
    "20D",
    "20E",
    "20G",
    "20H",
    "20I",
    "20J",
    "21A",
    "21B",
    "21C",
    "21D",
    "21E",
    "21F",
    "21G",
    "21H",
    "21I",
    "21J",
    "21K",
    "21L",
    "21M",
    "22A",
    "22B",
    "22C",
    "22D",
    "22E",
    "22F",
    "23A"
  ),
  `VOC/VOI` = c(
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Pre-VOC",
    "Alpha",
    "Gamma",
    "Delta",
    "Kappa",
    "Epsilon",
    "Eta",
    "Theta",
    "Iota",
    "Lambda",
    "Mu",
    "Delta",
    "Delta",
    "Omicron BA.1",
    "Omicron BA.2",
    "Omicron other",
    "Omicron BA.4",
    "Omicron BA.5",
    "Omicron BA.2.12.1",
    "Omicron BA.2.75",
    "Omicron BQ.1.1",
    "Omicron XBB.1",
    "Omicron XBB.1"
  )
)

data_kiemsurv <- data_kiemsurv_org %>%
  filter(Monsterstroom %in% c("TESTSTRAAT")) %>%
  left_join(relational_table_clade_VOC,
            by = "Clade") %>%
  mutate(
    `Test date` = as.Date(`Datum-monstername`, format = "%d-%m-%Y"),
    `Test week` = floor_date(`Test date`, unit = "week", week_start = 1),
    VOC = if_else(
      `VOC/VOI` %in% c(
        "Pre-VOC",
        "Alpha",
        "Beta",
        "Gamma",
        "Delta",
        "Omicron BA.1",
        "Omicron BA.2",
        "Omicron BA.4",
        "Omicron BA.5",
        "Omicron BA.2.12.1",
        "Omicron BA.2.75",
        "Omicron BQ.1.1",
        "Omicron XBB.1"
      ),
      `VOC/VOI`,
      NA_character_
    ) %>%
      factor(
        levels = c(
          "Pre-VOC",
          "Alpha",
          "Beta",
          "Gamma",
          "Delta",
          "Omicron BA.1",
          "Omicron BA.2",
          "Omicron BA.4",
          "Omicron BA.5",
          "Omicron BA.2.12.1",
          "Omicron BA.2.75",
          "Omicron BQ.1.1",
          "Omicron XBB.1"
        )
      )
  )

koppeltabel_variant_S_result <- bind_rows( # nolint
  expand.grid(
    `S result` = "Detected",
    variant_sgtf = "Pre-VOC",
    date = seq(
      from = as_date("2021-01-18"),
      to = as_date("2021-02-17"),
      by = "1 day"
    )
  ),
  expand.grid(
    `S result` = "Not detected",
    variant_sgtf = "Alpha",
    date = seq(
      from = as_date("2021-01-18"),
      to = as_date("2021-09-27"),
      by = "1 day"
    )
  ),
  expand.grid(
    `S result` = "Detected",
    variant_sgtf = "Delta",
    date = seq(
      from = as_date("2021-06-20"),
      to = as_date("2022-01-07"),
      by = "1 day"
    )
  ),
  expand.grid(
    `S result` = "Not detected",
    variant_sgtf = "Omicron BA.1",
    date = seq(
      from = as_date("2021-11-23"),
      to = as_date("2022-04-09"),
      by = "1 day"
    )
  ),
  expand.grid(
    `S result` = "Detected",
    variant_sgtf = "Omicron BA.2",
    date = seq(
      from = as_date("2022-01-29"),
      to = as_date("2022-10-18"),
      by = "1 day"
    )
  ),
  expand.grid(
    `S result` = "Not detected",
    variant_sgtf = "Omicron BA.4/5",
    date = seq(
      from = as_date("2022-05-02"),
      to = as_date("2022-10-05"),
      by = "1 day"
    )
  )
)

# tabel_ppv <- data_kiemsurv_sgtf %>%
#   filter(!is.na(VOC) & `S result` %in% c("Detected", "Not detected")) %>%
#   count(`S result`, VOC, `Datum-monstername`) %>%
#   group_by(`S result`, `Datum-monstername`) %>% 
#   mutate(
#     ppv = n / sum(n) * 100
#   ) %>% 
#   arrange(`Datum-monstername`, `S result`) %>% view

data_kiemsurv <- data_kiemsurv %>%
  slice(1, .by = `CoronIT Id`)

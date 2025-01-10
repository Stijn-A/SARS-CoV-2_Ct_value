# import data

data_sgtf_org <-
  read_rds(
    "/PATH/SGTF_data_20230221_1334.rds" # nolint
  )

source(
  "/PATH/functie_laad_data.R" # nolint
)

data_teststr_org <-
  functie_laad_data(
    "/PATH/Data/Teststraten/Geschoond/",
    as_date("2023-03-21"),
    fstcolumns = c(
      "Monsternummer",
      "Pseudoniem",
      "Afspraak_start_datum",
      "Uitslag",
      "Leeftijd",
      "Geslacht",
      "Is_gevaccineerd",
      "Vaccinatie_status",
      "Vaccinatie_datum_laatste",
      "Vaccinatie_aantal",
      "Vaccinatie_merk",
      "Vaccinatie_DUFS_EZD",
      "Ziektegegevens_eerste_ziektedag",
      "Klacht_Minstens_een"
    )
  ) %>%
  rename(
    Age = Leeftijd,
    sex = Geslacht,
    Is_vaccinated = Is_gevaccineerd,
    Vaccination_status = Vaccinatie_status,
    Vaccination_date_latest = Vaccinatie_datum_laatste,
    Vaccination_nr_dose = Vaccinatie_aantal,
    Vaccination_brand = Vaccinatie_merk,
    DUFS_EZD_0_42 = Vaccinatie_DUFS_EZD,
    Onset_date = Ziektegegevens_eerste_ziektedag,
    Symptoms_at_least_1 = Klacht_Minstens_een
  )

data_kiemsurv_org <-
  "/PATH/Covid_variants/data/" %>%
  list.files(full.names = TRUE) %>%
  str_subset("Lijst_voor_Epi") %>%
  sort() %>%
  str_subset("20230216.xlsx") %>%
  read_xlsx()

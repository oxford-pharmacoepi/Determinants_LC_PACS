# Load data ----
baselineCharacteristics <- read.csv(paste0(dir_results,"1-LoadCovariates/baselineCharacteristics.csv"))
biomarkers <- read.csv(paste0(dir_results,"1-LoadCovariates/biomarkers.csv"))
comorbidities <- read.csv(paste0(dir_results,"1-LoadCovariates/commorbidities.csv"))

# Baselines tables -----
baselineCharacteristics <- baselineCharacteristics |>
  addCoding("sex") |>
  addCoding("ethnic_background") |>
  addCoding("smoking_status")

baselineCharacteristics <- baselineCharacteristics |>
  summarise("N_total"  = prettyNum(paste0(nrow(baselineCharacteristics)), big.mark = ","),

            "Sex (%)_sex"  = "",
            "\tFemale_sex" = bin(sex, "Female"),
            "\tMale_sex"   = bin(sex, "Male"),
            "Counts_sex"   = counts(sex),
            "Missings_sex" = miss(sex),

            "Year of birth_year of birth"  = "",
            "\tMean (SD)_year of birth" = cont(year_of_birth),
            "\tQuantiles (0.05, 0.25, 0.5, 0.75, 0.95)_year of birth" = quantiles(year_of_birth),
            "Counts_year of birth"   = counts(year_of_birth),
            "Missings_year of birth" = miss(year_of_birth),

            "BMI (kg/m2)_body mass index"   = "",
            "\tMean (SD)_body mass index" = cont(body_mass_index),
            "\tQuantiles (0.05, 0.25, 0.5, 0.75, 0.95)_body mass index" = quantiles(body_mass_index),
            "Counts_body mass index" =  counts(body_mass_index),
            "Missings_body mass index" = miss(body_mass_index),

            "Index of multiple deprivation_index of multiple deprivation" = "",
            "\tMean (SD)_index of multiple deprivation" = cont(index_of_multiple_deprivation),
            "\tQuantiles (0.05, 0.25, 0.5, 0.75, 0.95)_index of multiple deprivation" = quantiles(index_of_multiple_deprivation),
            "Counts [N (%)]_index of multiple deprivation" =  counts(index_of_multiple_deprivation),
            "Missings_index of multiple deprivation" = miss(index_of_multiple_deprivation),

            "Ethnic background (%)_ethnic background" = "",
            "\tWhite_ethnic background" = bin(ethnic_background,"white"),
            "\t\tBritish_ethnic background" = bin(ethnic_background, "British"),
            "\t\tIrish_ethnic background" = bin(ethnic_background, "Irish"),
            "\t\tAny other white background_ethnic background" = bin(ethnic_background, "Any other white background"),

            "\tMixed_ethnic background_ethnic background" = bin(ethnic_background,"Mixed"),
            "\t\tWhite and Black Caribbean_ethnic background" = bin(ethnic_background, "White and Black Caribbean"),
            "\t\tWhite and Black African_ethnic background" = bin(ethnic_background, "White and Black African"),
            "\t\tWhite and Asian_ethnic background" = bin(ethnic_background, "White and Asian"),
            "\t\tAny other mixed background_ethnic background" = bin(ethnic_background, "Any other mixed background"),

            "\tAsian or Asian British_ethnic background" = bin(ethnic_background,"Asian or Asian British"),
            "\t\tIndian_ethnic background" = bin(ethnic_background, "Indian"),
            "\t\tPakistani_ethnic background" = bin(ethnic_background, "Pakistani"),
            "\t\tBangladeshi_ethnic background" = bin(ethnic_background, "Bangladeshi"),
            "\t\tAny other Asian background_ethnic background" = bin(ethnic_background, "Any other Asian background"),

            "\tBlack or Black British_ethnic background_ethnic background" = bin(ethnic_background,"Black or Black British"),
            "\t\tCaribbean_ethnic background" = bin(ethnic_background, "Caribbean"),
            "\t\tAfrican_ethnic background" = bin(ethnic_background, "African"),
            "\t\tAny other Black background_ethnic background" = bin(ethnic_background, "Any other Black background"),

            "\tChinese_ethnic background" = bin(ethnic_background, "Chinese"),
            "\tOther ethnic group_ethnic background" = bin(ethnic_background, "Other ethnic group"),

            "Counts_ethnic background_ethnic background" = counts(ethnic_background),
            "Missings_ethnic background_ethnic background" = miss(ethnic_background),

            "Smoking status (%)_smoking status" = "",
            "\tNever_smoking status" = bin(smoking_status, "Never"),
            "\tPrevious_smoking status" = bin(smoking_status, "Previous"),
            "\tCurrent_smoking status" = bin(smoking_status, "Current"),
            "Counts_smoking status" = counts(smoking_status),
            "Missings_smoking status" = miss(smoking_status)) |>
  mutate_all(~as.character(.)) |>
  pivot_longer(everything(), names_to = "Risk factor", values_to = "UK Biobank dataset") |>
  mutate(Name = gsub(" ","_",gsub(".*_","",`Risk factor`))) |>
  mutate(`Risk factor` = gsub("_.*","",`Risk factor`)) |>
  rowwise() |>
  mutate("z" = if_else(
    Name == "total",
    list(baselineCharacteristics[["sex"]][!is.na(baselineCharacteristics[["sex"]])]),
    list(baselineCharacteristics[[Name]][!is.na(baselineCharacteristics[[Name]])]))) |>
  ungroup() |>
  group_by(Name) |>
  mutate(z = if_else(row_number() == 1, z, list(NULL))) |>
  ungroup()

# Biomarkers -----
biomarkers <- summariseBiomarkers(biomarkers)

# Comorbidities ----
comorbidities <- comorbidities |>
  loadComorbidities() |>
  summariseComorbidities()

# Merge tables ----
baselineCharacteristics |>
  select("Risk factor", "UK Biobank dataset", "Name") |>
  filter(`Risk factor` == "N") |>
  rbind(
    tibble(
      "Risk factor" = "Baseline characteristics",
      "UK Biobank dataset" = " ",
      "Name" = "baseline_characteristics"
    )) |>
  rbind(
    baselineCharacteristics |> select("Risk factor", "UK Biobank dataset", "Name") |> filter(`Risk factor` != "N")
  ) |>
  rbind(
    tibble(
      "Risk factor" = "Biomarkers",
      "UK Biobank dataset" = " ",
      "Name" = "biomarkers"
    )
  ) |>
  rbind(
    biomarkers |> select("Risk factor", "UK Biobank dataset", "Name")
  ) |>
  rbind(
    tibble(
      "Risk factor" = "Comorbidities",
      "UK Biobank dataset" = " ",
      "Name" = "comorbidities"
    )
  ) |>
  rbind(
    comorbidities |> select("Risk factor", "UK Biobank dataset", "Name")
  ) |>
  rbind(
    tibble(
      "Risk factor" = "Vaccination status",
      "UK Biobank dataset" = " ",
      "Name" = "vaccine"
    )
  ) |>
  createFlextable(name = "ExploratoryAnalysis_original")


# Get figures ----
baselineCharacteristics |>
  mutate(type = if_else(`Risk factor` %in% c("N", "Sex (%)", "Ethnic background (%)", "Smoking status (%)"), "histogram", "line")) |>
  getFigures("Original")

baselineCharacteristics |>
  mutate(type = if_else(`Risk factor` %in% c("N", "Sex (%)", "Ethnic background (%)", "Smoking status (%)"), "histogram", "boxplot")) |>
  getFigures("Original")

biomarkers |> mutate(type = "line") |> getFigures("Original")

comorbidities |> mutate(type = "histogram") |> getFigures("Original")

rm(list = c("baselineCharacteristics", "biomarkers", "cohort", "comorbidities"))

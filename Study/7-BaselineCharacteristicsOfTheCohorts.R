# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |> select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |> select(-c("X"))

# Table ones ----

# Long covid
x <- read_csv(paste0(dir_data, "Results/5-LongCovid_cohort.csv"))
name <- "Long COVID"
name_cohort <- c("Controls","Cases")
longCovid_cohort <- tableOneStep1(x, baselineCharacteristics, biomarkers, name, name_cohort)

# PACS
x <- read_csv(paste0(dir_data, "Results/6-Pacs_cohort.csv"))
name <- "PACS"
name_cohort <- c("Controls","Cases")
pacs_cohort <- tableOneStep1(x, baselineCharacteristics, biomarkers, name, name_cohort)

# Merge both tables ----
y <- longCovid_cohort |>
  mutate(order = row_number()) |>
  inner_join(pacs_cohort |>  mutate(order = row_number()), by = c("Risk factor", "order")) |>
  select(-c("order")) |>
  mutate(`Risk factor` = if_else(`Risk factor` == "Sociodemographics factors", "Sociodemographic factors",`Risk factor`))

y <- y |>
  filter(`Risk factor` == "N") |>
  rbind(
    y |>
      filter(`Risk factor` != "N")
  )

merge <- c(1:nrow(y))[y$`Risk factor` %in% c("Sociodemographic factors", "Comorbidities [Cases (%)]", "Biomarkers [Mean (SD)]")]
row_fill <- c(1:nrow(y))[!y$`Long COVID_Controls` == " "]
y <- y |>
  flextable() |>
  span_header() |>
  bold(i = 1, part = "header") |>
  bold(i = 2, part = "header") |>
  align(j = c(2,3,4,5), align = "center", part = "all") |>
  width(j = 1, width = 5, unit = "cm") |>
  bg(bg = "#F2F2F2", i = 1, j = c(2,3), part = "header") |>
  bg(bg = "#F2F2F2", i = 2, j = c(2,4), part = "header") |>
  bg(bg = "#F2F2F2", i = row_fill, j = c(2,4), part = "body")

for(ii in merge){
  y <- y |>
    merge_at(i = ii, j = c(1:5)) |>
    hline(i = ii-1, border = fp_border(color = "black")) |>
    bold(i = ii)
}

y |>
  save_as_docx(path = paste0(dir_results, "Tables/BaselineCharacteristics.docx"))

rm(list = c("baselineCharacteristics", "biomarkers","x","name","name_cohort",
            "longCovid_cohort","pacs_cohort","y", "merge", "row_fill"))

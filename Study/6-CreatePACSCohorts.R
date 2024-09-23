# Load databases ---------------------------------------------------------------
covid19_result <- loadCovid19Result()
hes_data       <- loadHesData()

# Load phenotypes --------------------------------------------------------------
sequela_table <- loadSequelaTable()

# Build cohort -----------------------------------------------------------------
t <- hes_data |>
  recordAttrition() |>
  mutate(diag_icd10   = if_else(!diag_icd10 %in% sequela_table$icd10_code, NA, diag_icd10),
         episode_date = if_else(is.na(diag_icd10), NA, episode_date)) |>
  distinct() |>
  recordAttrition("Restrict to PACS records") |>
  inner_join(
    covid19_result |>
      select("eid","specdate","result"),
    by = "eid",
    relationship = "many-to-many"
  ) |>
  recordAttrition("Restrict to participants with COVID-19 linkage data") |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  mutate(episode_date = if_else(episode_date >= (specdate-365), episode_date, NA)) |>
  mutate(episode_date = if_else(episode_date <= (specdate+365), episode_date, NA)) |>
  mutate(diag_icd10 = if_else(is.na(episode_date), NA, diag_icd10)) |>
  distinct() |>
  recordAttrition("Restrict the study period between one year before and after testing positive for COVID-19.") |>
  mutate(exclusion = if_else(
    episode_date >= (specdate-365) & episode_date <= specdate & (!is.na(diag_icd10)),
    1,0
  ))

t1 <- t |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t) |>
  filter(exclusion == 0) |>
  recordAttrition("Restrict to participants that did not have any PACS event one year prior the COVID-19 infection") |>
  mutate(exclusion = if_else((!is.na(diag_icd10)) &
                               (episode_date > specdate) &
                               (episode_date <= (specdate+washout_period)), 1, 0))

t2 <- t1 |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t1) |>
  filter(exclusion == 0) |>
  recordAttrition(paste0("Restrict to participants that did not have any PACS event after ", washout_period, " days of the infection")) |>
  mutate(diagnoses = if_else(!is.na(diag_icd10), 1, 0))

t3 <- t2 |>
  group_by(eid) |>
  mutate(cases = sum(diagnoses)) |>
  ungroup() |>
  restoreAttrition(t2)

pacs_cases_cohort <- t3 |>
  filter(cases >= n_symptoms) |>
  filter(!is.na(diag_icd10)) |>
  recordAttrition("Restrict to cases")

pacs_controls_cohort <- t3 |>
  filter(cases >= 0 & cases < n_symptoms) |>
  distinct() |>
  recordAttrition("Restrict to controls")

pacs_cohort  <- pacs_cases_cohort |>
  union_all(pacs_controls_cohort) |>
  select("eid", "specdate", "state" = "cases") |>
  group_by(eid) |>
  mutate(specdate = min(specdate),
         state = if_else(state < n_symptoms, 0, 1)) |>
  distinct() |>
  ungroup()




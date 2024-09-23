# Load databases ----
health_questionnaire <- loadHealthAndWellBeingQuestionnaire()
covid19_result       <- loadCovid19Result()

# Create cohort ----
x <- covid19_result |>
  select("eid","specdate","result") |>
  recordAttrition() |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  inner_join(
    health_questionnaire |>
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) |>
  recordAttrition("Restrict to participants that answered the health and well-being web-questionnaire") |>
  filter(questionnaire_started != as.Date("1999-01-01")) |>
  recordAttrition("Restrict to participants that answered yes/no to all the questions from the health and well-being web-questionnaire.") |>
  filter(specdate < questionnaire_started) |>
  recordAttrition("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction before answering the questionnaire.")

x1 <- x |> group_by(eid) |>
  filter(specdate == max(specdate)) |>
  ungroup() |>
  restoreAttrition(x) |>
  recordAttrition("Restrict to the latest record.")

x2 <- x1 |>
  filter(((specdate+washout_period) < questionnaire_started) & (specdate > (questionnaire_started-365))) |>
  recordAttrition(paste0("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction between 1 year and ",washout_period," days before answering the questionnaire."))

x3 <- x2 |>
  rowwise() |>
  mutate("symptom" = across(starts_with("symptom")) |> sum()) |>
  ungroup() |>
  restoreAttrition(x2) %>%
  mutate(length = do.call(pmax, c(select(., starts_with("length")), na.rm = TRUE))) |>
  select(-starts_with("length_")) |>
  filter((questionnaire_started-length) > specdate) |>
  recordAttrition("Restrict to people with no persistent symptoms")

longCovid_cases_cohort <- x3 |>
  filter(symptom >= n_symptoms) |>
  recordAttrition("Restrict to people that reported at least one symptom of the questionnaire")

longCovid_controls_cohort <- x3 |>
  filter(symptom < n_symptoms & symptom >= 0) |>
  recordAttrition("Restrict to people that did not report any symptom of the questionnaire")

longCovid_cohort <- longCovid_cases_cohort |>
  select("eid", "specdate", "questionnaire_started") |>
  mutate(state = 1) |>
  union_all(
    longCovid_controls_cohort |>
      select("eid", "specdate", "questionnaire_started") |>
      mutate(state = 0)
  )







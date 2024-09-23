# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |> select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |> select(-c("X"))

nonLinear <- as_tibble(read.csv(paste0(dir_results,"9-LinearityAnalysis/", cohort_name, "_list.csv"))) |>
  select("name", "significance") |> distinct() |>
  filter(significance == 1)

# Merge datasets ----
cohort <- mergeAllCovariates(read_csv(paste0(dir_data, "Results/",cohort_name,"_cohort.csv")),
                             baselineCharacteristics,
                             biomarkers)

if(analysis == "secondary"){
  state <- as_tibble(read.csv(paste0(dir_results,"SecondaryAnalysis/",cohort_name,"_cohort.csv")))
  cohort <- cohort |> select(-"state") |> inner_join(state |> select("eid","state"), by = "eid")
}

# Lasso features ----
lassoFeatures <- list()
lassoFeatures[[1]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/baselineCharacteristics_",cohort_name,".csv"))
lassoFeatures[[2]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/biomarkers_",cohort_name,".csv"))
lassoFeatures[[3]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/comorbidities_",cohort_name,".csv"))

# Data to adjust by ----
data <- list()
data[[1]] <- cohort |>
  select(lassoFeatures[[1]]$Lasso.features,"state", "eid") |>
  reverseCoding("sex") |>
  reverseCoding("smoking_status") |>
  mutate(ethnic_background = if_else(ethnic_background == "White", 0, 1))

data[[2]] <- cohort |>
  select(gsub("--.*","",lassoFeatures[[2]]$Lasso.features) |> unique(), "state", "eid") |>
  select(-c("sex", "age_when_infected"))

data[[3]] <- cohort |>
  select(lassoFeatures[[3]]$Lasso.features, "state", "eid") |>
  select(-c("age_when_infected","sex")) |>
  mutate(across(-c("state", "eid"), ~if_else(. == "Cases", 1,0)))

tableResults <-tibble("Estimate"   = as.double(),
                      "Std. Error" = as.double(),
                      "z value"    = as.double(),
                      "Pr(>|z|)"   = as.double(),
                      "Risk factor" = as.character(),
                      "Variable" = as.character(),
                      "Type" = as.character())
results <- list()

# Categorise nonlinear biomarkers ----
for(var in  colnames(data[[2]])[!colnames(data[[2]]) %in% c("state", "age_when_infected", "sex")]){
  x1 <- data[[2]] |> select("risk_factor" = var)

  if(var %in% nonLinear$name){
    x1 <- categoriseBiomarker(x1)
  }

  data[[2]][[var]] <- x1$risk_factor
}


# Merge all data ----
dataAnalysis <- data[[1]] |>
  categoriseVariable("body_mass_index") |>
  categoriseVariable("index_of_multiple_deprivation") |>
  categoriseVariable("age_when_infected") |>
  mutate(smoking_status = as.factor(smoking_status)) |>
  inner_join(
    data[[3]], by = c("eid","state")
  ) |>
  inner_join(
    data[[2]], by = c("eid","state")
  ) |>
  select(-c("eid","age"))


# Adjusted analysis
result <- glm(data = dataAnalysis, family = "binomial", formula = state ~ .)

# Check multicollinearity
vif_before <- vif(result) |>
  as_tibble() |>
  mutate("Risk factor" = vif(result) |> row.names()) |>
  select("Risk factor", "GVIF") |>
  left_join(
    tibbleNameComorbidities() |>
      full_join(tibbleNameBiomarkers()) |>
      rename("Risk factor"  = "risk_factor",
             "Risk factor1" = "risk_factor1"),
    by = "Risk factor") |>
  mutate(`Risk factor1` = if_else(is.na(`Risk factor1`), str_to_sentence(gsub("_"," ",`Risk factor`)), `Risk factor1`)) |>
  select("Risk factor" = "Risk factor1", "GVIF") |>
  mutate(GVIF = round(GVIF,2))

result <- glm(data = dataAnalysis, family = "binomial", formula = state ~ .)
summary(result)$coefficients |>
  as_tibble() |>
  mutate(Name = row.names(summary(result)$coefficients)) |>
  filter(Name != "(Intercept)") |>
  mutate(n = ifelse(substr(Name, nchar(Name), nchar(Name)) %in% c(1,2,3,4,5) & Name != "IGF1",
                    substr(Name, nchar(Name), nchar(Name)),
                    "1")) |>
  mutate(Name = ifelse(substr(Name, nchar(Name), nchar(Name)) %in% c(1,2,3,4,5)  & Name != "IGF1",
                       substr(Name, 1, nchar(Name)-1),
                       Name)) |>
  full_join(
    nameCategorised() |>
      rename("Name" = `Risk factor`) |> mutate(n = as.character(n)),
    by = c("Name", "n")
  ) |>
  group_by(Name) |>
  mutate(Variable = if_else(is.na(Variable) & n() > 1,
                            paste0(" - Q",n),
                            Variable)) |>
  mutate(OR = exp(Estimate),
         CI_LOW = exp(Estimate-1.96*`Std. Error`),
         CI_HIGH = exp(Estimate+1.96*`Std. Error`)) |>
  select(`Risk factor` = "Name", Variable, OR, CI_LOW, CI_HIGH, `Pr(>|z|)`) |>
  mutate(Analysis = "Adjusted") |>
  write.csv(if_else(analysis == "secondary",
                    paste0(dir_results,"SecondaryAnalysis/",cohort_name, "_LogisticResults_Adjusted.csv"),
                    paste0(dir_results,cohort_name, "_LogisticResults_Adjusted.csv")))

vif(result) |>
  as_tibble() |>
  mutate("Risk factor" = vif(result) |> row.names()) |>
  select("Risk factor", "GVIF") |>
  left_join(
    tibbleNameComorbidities() |>
      full_join(tibbleNameBiomarkers() |> arrange(risk_factor1),
                by = c("risk_factor","risk_factor1")) |>
      rename("Risk factor" = "risk_factor", "Risk factor1" = "risk_factor1"),
    by = "Risk factor") |>
  mutate(`Risk factor1` = if_else(is.na(`Risk factor1`), str_to_sentence(gsub("_"," ",`Risk factor`)), `Risk factor1`)) |>
  select("Risk factor" = "Risk factor1", "GVIF_after" = "GVIF") |>
  mutate(GVIF_after = round(GVIF_after,2)) |>
  select(`Risk factor`,"GVIF (after)" = "GVIF_after") |>
  flextable() |>
  bold(part = "header") |>
  width(j = 1, width = 6, unit = "cm") |>
  save_as_docx(path = paste0(dir_results, "Tables/", cohort_name, "VarianceInflationFactor.docx"))


rm(list = c("baselineCharacteristics","biomarkers","nonLinear","cohort","lassoFeatures","data",
            "tableResults","results","baselineCharacteristics_results","names","var","x1","results",
            "tableResults","comorbidities_results","biomarkers_results"))



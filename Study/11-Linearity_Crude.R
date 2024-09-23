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

cohort <- cohort |> select(-c("eid"))

# Lasso features ----
lassoFeatures <- list()
lassoFeatures[[1]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/baselineCharacteristics_",cohort_name,".csv"))
lassoFeatures[[2]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/biomarkers_",cohort_name,".csv"))
lassoFeatures[[3]] <- read.csv(paste0(dir_data,"Results/10-LassoRegression/comorbidities_",cohort_name,".csv"))

# Data to adjust by ----
data <- list()
data[[1]] <- cohort |>
  select(lassoFeatures[[1]]$Lasso.features,"state") |>
  reverseCoding("sex") |>
  reverseCoding("smoking_status") |>
  mutate(ethnic_background = if_else(ethnic_background == "White", 0, 1))

data[[2]] <- cohort |>
  select(gsub("--.*","",lassoFeatures[[2]]$Lasso.features) |> unique(), "state") |>
  reverseCoding("sex")

data[[3]] <- cohort |>
  select(lassoFeatures[[3]]$Lasso.features, "state") |>
  mutate(across(-c("age_when_infected","sex","state"), ~if_else(. == "Cases", 1,0))) |>
  reverseCoding("sex")

tableResults <-tibble("Estimate"   = as.double(),
                      "Std. Error" = as.double(),
                      "z value"    = as.double(),
                      "Pr(>|z|)"   = as.double(),
                      "Risk factor" = as.character(),
                      "Variable" = as.character(),
                      "Type" = as.character())
results <- list()

# BaselineCharacteristics results ----
baselineCharacteristics_results <- data[[1]] |>
  categoriseVariable("body_mass_index") |>
  categoriseVariable("index_of_multiple_deprivation") |>
  categoriseVariable("age_when_infected") |>
  mutate(smoking_status = as.factor(smoking_status))

names <- colnames(baselineCharacteristics_results)[!colnames(baselineCharacteristics_results) %in% c("state", "age")]

for(var in names){
  x1 <- baselineCharacteristics_results |> rename("risk_factor" = var)

  if(var == "sex"){
    results[[var]] <- glm(data = x1, state ~ age + risk_factor, family = "binomial")
  }else if(var == "age_when_infected"){
    results[[var]] <- glm(data = x1, state ~ sex + risk_factor, family = "binomial")
  }else{
    results[[var]] <- glm(data = x1, state ~ age + sex + risk_factor, family = "binomial")
  }

  tableResults <- tableResults |>
    rbind(
      summary(results[[var]])$coefficients |>
        as_tibble() |>
        mutate(Name = row.names(summary(results[[var]])$coefficients)) |>
        mutate("Risk factor" = c(var)) |>
        filter(row_number() > if_else(var %in% c("sex", "age_when_infected"), 2, 3)) |>
        mutate("n" = as.numeric(gsub("risk_factor","",Name))) |>
        mutate("n" = if_else(is.na(n), 1, n)) |>
        inner_join(
          nameCategorised(),
          by = c("Risk factor", "n")
        ) |>
        select(-c("n","Name")) |>
        mutate(Type = "Baseline characteristics")
    )
}

# Comorbidities results ----
comorbidities_results <- data[[3]]

names <- colnames(comorbidities_results)[!colnames(comorbidities_results) %in% c("state", "age_when_infected", "sex")]

for(var in names){
  x1 <- comorbidities_results |> rename("risk_factor" = var)

  results[[var]] <- glm(data = x1, state ~ age_when_infected + sex + risk_factor, family = "binomial")

  tableResults <- tableResults |>
    rbind(
      summary(results[[var]])$coefficients |>
        as_tibble() |>
        mutate("Risk factor" = c(var), "Variable" = "-") |>
        filter(row_number() == 4) |>
        mutate(Type = "Comorbidities")
    )
}

# Biomarkers results ----
biomarkers_results <- data[[2]]

names <- colnames(biomarkers_results)[!colnames(biomarkers_results) %in% c("state", "age_when_infected", "sex")]

for(var in names){
  x1 <- biomarkers_results |> rename("risk_factor" = var)

  if(var %in% nonLinear$name){
    x1 <- categoriseBiomarker(x1)
  }

  results[[var]] <- glm(data = x1, state ~ age_when_infected + sex + risk_factor, family = "binomial")

  tableResults <- tableResults |>
    rbind(
      summary(results[[var]])$coefficients |>
        as_tibble() |>
        mutate("Risk factor" = c(var)) |>
        filter(row_number() > 3) |>
        group_by(`Risk factor`) |>
        mutate(n = n()) |>
        ungroup() |>
        mutate("Variable" = if_else(n > 1, paste0("Q", row_number()+1), "")) |>
        select(-"n") |>
        mutate(Type = "Biomarkers")
    )
}


tableResults |>
  mutate(order = row_number()) |>
  arrange(desc(order)) |>
  mutate(order = row_number()) |>
  mutate(OR = exp(Estimate),
         CI_LOW = exp(Estimate-1.96*`Std. Error`),
         CI_HIGH = exp(Estimate+1.96*`Std. Error`)) |>
  select(`Risk factor`, Variable, OR, CI_LOW, CI_HIGH, `Pr(>|z|)`, Type) |>
  mutate(Analysis = "Crude") |>
  write.csv(if_else(analysis == "secondary",
                    paste0(dir_results,"SecondaryAnalysis/",cohort_name, "_LogisticResults_Crude.csv"),
                    paste0(dir_results,cohort_name, "_LogisticResults_Crude.csv")))

rm(list = c("baselineCharacteristics","biomarkers","nonLinear","cohort","lassoFeatures","data",
            "tableResults","results","baselineCharacteristics_results","names","var","x1","results",
            "tableResults","comorbidities_results","biomarkers_results"))




# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |> select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |> select(-c("X"))

nonLinear <- as_tibble(read.csv(paste0(dir_results,"9-LinearityAnalysis/", cohort_name, "_list.csv"))) |>
  select("name", "significance") |> distinct() |>
  filter(significance == 1)

# Merge datasets ----
cohortData <- read_csv(paste0(dir_data, "Results/",cohort_name,"_cohort.csv"))
cohort <- mergeAllCovariates(cohortData, baselineCharacteristics, biomarkers) |> select(-c("eid"))

dataLasso <- list()
dataLasso[[1]] <- cohort |>
  select(any_of(colnames(baselineCharacteristics)),"age_when_infected", "state") |>
  reverseCoding("sex") |>
  reverseCoding("smoking_status") |>
  mutate(ethnic_background = if_else(ethnic_background == "White", 0, 1))

dataLasso[[2]] <- cohort |>
  select(any_of(colnames(biomarkers)), "state", "age_when_infected","sex") |>
  reverseCoding("sex") |>
  select(-c("apolipoprotein_b", "cholesterol", "apolipoprotein_a", "direct_bilirubin",
            "aspartate_aminotransferase", "gamma_glutamyltransferase", "creatinine", "glucose",
            "total_protein","albumin", "testosterone"))

dataLasso[[3]] <- cohort |> select(age_when_infected, state, sex, (tibbleNameComorbidities() |> pull("risk_factor"))) |>
  mutate(across(-c("age_when_infected","sex","state"), ~if_else(. == "Cases", 1,0))) |>
  reverseCoding("sex") |>
  select(-c("asthma","liver_disease_mild"))

# Lasso regression ----
set.seed(43)
lambdas  <- list()
results  <- list()
for(i in 1:3){
  x <- dataLasso[[i]] |> select(-c("state")) |> as_tibble()
  y <- dataLasso[[i]] |> select("state") |> pull()

  if(i == 2){
    for(j in c(1:nrow(nonLinear))){
      x1 <- x |>
        mutate(n = row_number()) |>
        select("n", nonLinear$name[[j]]) |>
        arrange(.data[[nonLinear$name[[j]]]])

       x <- x |>
         select(-.data[[nonLinear$name[[j]]]]) |>
         mutate(n = row_number()) |>
         inner_join(
           ns(x1 |> select(-"n") |> pull(), 5, intercept = FALSE) |>
             as_tibble() |>
             rename_all(~paste0(nonLinear$name[[j]], "--",.)) |>
             mutate(n = x1$n) |>
             arrange(n),
           by = "n"
         ) |>
         select(-"n")
    }
  }

  x <- x |> as.matrix()

  # Cross validation to select the lambda
  lasso_reg <- cv.glmnet(
    type.measure = "deviance",
    x = x,
    y = y,
    lambda = 10^seq(2, -4, by = -.001),
    standardize = FALSE,
    nfolds = 10,
    family = "binomial",
    alpha = 1
  )

  # Extract coefficients
  coef.lasso_reg <- coef(lasso_reg, s = lasso_reg$lambda.min)
  selectedLassoFeatures <- names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
  selectedLassoFeatures <- selectedLassoFeatures[selectedLassoFeatures != "(Intercept)"]
  selectedLassoCoefficients <- as.numeric(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])[-1]

  results[[i]] <- tibble(
    "Lasso features" = selectedLassoFeatures,
    "Lasso coefficients" = selectedLassoCoefficients
  )

  lambdas[[i]] <- lasso_reg$lambda.min
}

# Merge all the results ----
tibble("Lasso features" = "Baseline characteristics", "Lasso coefficients" = 1000) |>
  rbind(
    tibble("Lasso features" = c("age_when_infected", "sex")) |>
      right_join(results[[1]] |> arrange(`Lasso features`), by = "Lasso features")
    ) |>
  rbind(tibble("Lasso features" = "Biomarkers", "Lasso coefficients" = 1000)) |>
  rbind(
    tibble("Lasso features" = c("age_when_infected", "sex")) |>
      right_join(results[[2]] |> arrange(`Lasso features`), by = "Lasso features")
  ) |>
  rbind(tibble("Lasso features" = "Comorbidities", "Lasso coefficients" = 1000)) |>
  rbind(
    rbind(
      tibble("Lasso features" = c("age_when_infected", "sex")) |>
        right_join(results[[3]] |> arrange(`Lasso features`), by = "Lasso features")
    )
  ) |>
  mutate(`Lasso features` = str_to_sentence(gsub("_", " ", `Lasso features`))) |>
  mutate(`Lasso coefficients` = round(`Lasso coefficients`, digits = 4)) |>
  flextable() |>
  bold(~`Lasso coefficients` == 1000) |>
  bold(i = 1, part = "header") |>
  width(j = 1, width = 7, unit = "cm") |>
  width(j = 2, width = 4, unit = "cm") |>
  merge_at(i = 1, j = c(1,2), part = "body") |>
  merge_at(i = nrow(results[[1]])+2, part = "body") |>
  merge_at(i = nrow(results[[2]])+nrow(results[[1]])+3, part = "body") |>
  italic(i = 1, part = "body") |>
  italic(i = nrow(results[[1]])+2, part = "body") |>
  italic(i = nrow(results[[2]])+nrow(results[[1]])+3, part = "body") |>
  hline(i = nrow(results[[1]])+1, part = "body") |>
  hline(i = nrow(results[[2]])+nrow(results[[1]])+2, part = "body") |>
  save_as_docx(path = paste0(dir_results,"Tables/lasso_features_",cohort_name,".docx"))

tibble("Analysis" = c("Baseline characteristics", "Biomarkers", "Comorbidities"),
       "Lambdas" = Reduce(union_all,lambdas)) |>
  write.csv(paste0(dir_results,"10-LassoRegression/lambdas_",cohort_name,".csv"), row.names = FALSE)


results[[1]] |> write.csv(paste0(dir_results,"10-LassoRegression/baselineCharacteristics_",cohort_name,".csv"), row.names = FALSE)
results[[2]] |> write.csv(paste0(dir_results,"10-LassoRegression/biomarkers_",cohort_name,".csv"), row.names = FALSE)
results[[3]] |> write.csv(paste0(dir_results,"10-LassoRegression/comorbidities_",cohort_name,".csv"), row.names = FALSE)

rm(list = c("baselineCharacteristics","biomarkers","cohort_name","cohortData","cohort",
            "dataLasso","lambdas","results","i","x","x1", "y","lasso_reg","coef.lasso_reg",
            "selectedLassoFeatures","selectedLassoCoefficients", "nonLinear"))


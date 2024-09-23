# Crude analysis ------
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |>
  select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |>
  select(-c("X"))
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
lassoFeatures <- read.csv(paste0(dir_data,"Results/10-LassoRegression/biomarkers_",cohort_name,".csv"))

# Data to adjust by ----
data <- cohort |>
  select(gsub("--.*","",lassoFeatures$Lasso.features) |> unique(), "state") |>
  reverseCoding("sex") |>
  select(any_of(c("sex", "age_when_infected", "testosterone", "shbg", "state")))

tableResults <-tibble("Estimate"   = as.double(),
                      "Std. Error" = as.double(),
                      "z value"    = as.double(),
                      "Pr(>|z|)"   = as.double(),
                      "Risk factor" = as.character(),
                      "Variable" = as.character(),
                      "Type" = as.character())
results <- list()

# Biomarkers results ----
biomarkers_results <- data
names <- colnames(biomarkers_results)[!colnames(biomarkers_results) %in% c("state", "age_when_infected", "sex")]

for(i in c(0,1)){
  biomarkers_results1 <- biomarkers_results |> filter(sex == i)
  for(var in names){
    x1 <- biomarkers_results1 |> rename("risk_factor" = var)

    if(var %in% nonLinear$name){
      x1 <- categoriseBiomarker(x1)
    }

    results[[var]] <- glm(data = x1, state ~ age_when_infected + risk_factor, family = "binomial")

    tableResults <- tableResults |>
      rbind(
        summary(results[[var]])$coefficients |>
          as_tibble() |>
          mutate("Risk factor" = c(var)) |>
          filter(row_number() > 2) |>
          group_by(`Risk factor`) |>
          mutate(n = n()) |>
          ungroup() |>
          mutate("Variable" = if_else(n > 1, paste0("Q", row_number()+1), "")) |>
          select(-"n") |>
          mutate(Type = "Biomarkers") |>
          mutate(Sex = if_else(i == 0, "Female", "Male"))
      )
  }
}

crude <- tableResults |>
  mutate(`Risk factor` = str_to_upper(`Risk factor`)) |>
  mutate(order = row_number()) |>
  arrange(desc(order)) |>
  mutate(order = row_number()) |>
  mutate(OR = exp(Estimate),
         CI_LOW = exp(Estimate-1.96*`Std. Error`),
         CI_HIGH = exp(Estimate+1.96*`Std. Error`)) |>
  select(`Risk factor`, Variable, OR, CI_LOW, CI_HIGH, `Pr(>|z|)`, Type, Sex) |>
  mutate(Analysis = "Crude") |>
  group_by(`Risk factor`) |>
  mutate(`Risk factor` = if_else(row_number() == 2, "", `Risk factor`)) |>
  mutate(`OR (95%CI)` = paste0(round(OR,2)," (", round(CI_LOW,2),", ", round(CI_HIGH,2),")"),
         `P-Value` = if_else(`Pr(>|z|)` < 0.05, "< 0.05", as.character(round(`Pr(>|z|)`,2)))) |>
  mutate("                                                                 " = "")

tm <- forest_theme(
                   base_family = "Calibri",
                   base_size = 10,
                   refline_lty = "solid",
                   ci_pch = c(15),
                   ci_col = "#8E2723",
                   ci_Theight = 0,
                   colhead =  list(fg_params = list(fontsize = 12, fontfamily = "Calibri")))

p <- forest(crude[,c(1,8,12,10,11)],
            est = list(crude$OR),
            lower = list(crude$CI_LOW),
            upper = list(crude$CI_HIGH),
            ci_column = 3,
            ref_line = 1,
            x_trans = "log",
            xlim = c(0.5,1.5),
            theme = tm,
            ticks_at = c(0.5, 0.75, 1, 1.5, 2),
            xlab = "OR") |>
  add_border(part = "header", row = 1, where = "bottom")

p_sc <- get_scale(plot = p, width_wanted = 6, height_wanted = 4, unit = "in")
ggplot2::ggsave(filename = paste0(dir_results,"Figures/Analysis3/",cohort_name,"_forestplot.png"),
                plot = p,
                dpi = 500,
                width = 6,
                height = 4,
                units = "in",
                scale = p_sc)

rm(list = c("baselineCharacteristics","biomarkers","nonLinear","cohort","lassoFeatures","data",
            "tableResults","results","baselineCharacteristics_results","names","var","x1","results",
            "tableResults","comorbidities_results","biomarkers_results"))





# Adjusted analysis ----
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
  select(-c("sex", "age_when_infected")) |>
  select(-any_of(c("alanine_aminotransferase", "alkaline_phosphatase", "aspartate aminotransferase",
                   "gamma_glutamyltransferase", "hbauc",  "glucose", "igf")))

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
adjusted <- summary(result)$coefficients |>
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
  mutate(Analysis = "Adjusted")


rm(list = c("baselineCharacteristics","biomarkers","nonLinear","cohort","lassoFeatures","data",
            "tableResults","results","baselineCharacteristics_results","names","var","x1","results",
            "tableResults","comorbidities_results","biomarkers_results"))


# Plot results ----
windowsFonts("Calibri" = windowsFont("Calibri"))

get_scale <- function(plot,
                      width_wanted,
                      height_wanted,
                      unit = "in"){
  h <- convertHeight(sum(plot$heights), unit, TRUE)
  w <- convertWidth(sum(plot$widths), unit, TRUE)
  max(c(w/width_wanted,  h/height_wanted))
}

adjusted <- tibbleNameComorbidities() |>
  mutate(Type = "Comorbidities") |>
  full_join(tibbleNameBiomarkers() |>
              mutate(Type = "Biomarkers"),
            by = c("risk_factor", "risk_factor1", "Type")) |>
  full_join(
    adjusted |>
      rename("risk_factor" = "Risk factor", "PValue" = "Pr(>|z|)"),
    by = "risk_factor"
  ) |>
  filter((!is.na(OR)) | (is.na(OR) & risk_factor == risk_factor1)) |>
  mutate(Variable = gsub(" - ","",Variable)) |>
  mutate(Variable = if_else(is.na(Variable), "", Variable)) |>
  mutate(Variable = gsub("Q4","Q5", Variable),
         Variable = gsub("Q3","Q4", Variable),
         Variable = gsub("Q2","Q3", Variable),
         Variable = gsub("Q1","Q2", Variable)) |>
  mutate(Type = if_else(is.na(Type), "Baseline characteristics", Type))  |>
  mutate(Analysis = "Adjusted")


l1 <- c(5,3,5)
col <- c("#5AAEB2", "#BC5856", "#F3CC45")
nam <- c("Baseline characteristics", "Biomarkers", "Comorbidities")
nam1 <- c("A) Sociodemographic factors", "B) Biomarkers", "C) Comorbidities")
plist <- list()

for(i in 1:3){
  t_adj   <- adjusted |>
    filter(Type == nam[i]) |>
    mutate(risk_factor1 = if_else(is.na(risk_factor1), str_to_sentence(gsub("_"," ",risk_factor)), risk_factor1)) |>
    filter(Analysis == "Adjusted" | is.na(Analysis)) |>
    rename("OR_adj" = "OR", "CI_LOW_adj" = "CI_LOW", "CI_HIGH_adj" = "CI_HIGH",
           "PValue_adj" = "PValue")

  if(cohort_name == "6-Pacs"){
    t_adj <- t_adj |> filter(risk_factor != "Liver")
  }

  dt <- t_adj |>
    select("risk_factor" = "risk_factor1", "Type", "Variable", ends_with("_adj")) |>
    mutate(n = rev(row_number())) |>
    mutate("OR"   = if_else(is.na(OR_adj),"", paste0(round(OR_adj,2)," (", round(CI_LOW_adj,2),", ", round(CI_HIGH_adj,2),")")),
           "Pval" = if_else(is.na(PValue_adj),"",
                            if_else(PValue_adj < 0.05, "< 0.05",
                                    paste0(formatC(PValue_adj, format = "e", 1)))))
  if(i == 1){
    dt <- dt |>
      group_by(risk_factor) |>
      mutate(n1 = row_number()) |>
      arrange(risk_factor, n1) |>
      select(-"n1") |>
      ungroup()
  }

  if(i == 2){  dt <- dt |> filter(risk_factor != "Diabetes")}

  if(i == 2){
    dt <- dt |>
      mutate("                                                                                                                                  " = "")
  }else{
    dt <- dt |>
      mutate("                                                                 " = "")
  }

  dt <- dt |>
    mutate(risk_factor = if_else(row_number() > 1, "",risk_factor), .by = c("risk_factor")) |>
    mutate(risk_factor = if_else(risk_factor == "Age when infected", "Age when infected\n(Ref: <55)",risk_factor)) |>
    mutate(risk_factor = if_else(risk_factor == "Body mass index", "Body mass index\n(Ref: Average)",risk_factor)) |>
    mutate(risk_factor = if_else(risk_factor == "Ethnic background", "Ethnic background\n(Ref: White)",risk_factor)) |>
    mutate(risk_factor = if_else(risk_factor == "Index of multiple deprivation", "Index of multiple deprivation\n(Ref: Very affluent)",risk_factor)) |>
    mutate(risk_factor = if_else(risk_factor == "Sex", "Sex\n(Ref: Women)",risk_factor)) |>
    mutate(risk_factor = if_else(risk_factor == "Smoking status", "Smoking status\n(Ref: Never)",risk_factor)) |>
    rename(!!nam1[i] := "risk_factor",
           " " = "Variable",
           "OR (95% CI)" = "OR",
           "P-Value" = "Pval")

  rm("t_adj")

  plist[[i]] <- getUngroupedForestPlot(dt, dir_results, cohort_name, nam[i], l1[i], col[i])

}

max_width <- unit.pmax(plist[[1]]$widths, plist[[3]]$widths, plist[["2.1"]]$widths, plist[["2.2"]]$widths)

# Set the same widths for both plots
i <- 1
for(p in plist){
  if(i == 2 && cohort_name == "6-Pacs"){
    h <- 15
    w <- 10
  }else{
    h <- 7
    w = 8
  }
  p$widths <- max_width
  p_sc <- get_scale(plot = p, width_wanted = w, height_wanted = h, unit = "in")
  ggplot2::ggsave(filename = paste0(dir_results,"Figures/Analysis2/",cohort_name,"_forestplot_main_",i,".png"),
                  plot = p,
                  dpi = 500,
                  units = "in",
                  width = w, height = h,
                  scale = p_sc
  )
  i <- i + 1
}


rm(list = c("crude","adjusted", "dataAnalysis","t","t1","i","l1","nam","path"))


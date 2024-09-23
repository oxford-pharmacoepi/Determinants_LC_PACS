# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |> select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |>
  select(-c("X")) |>
  select(-c("apolipoprotein_b", "cholesterol", "apolipoprotein_a", "direct_bilirubin",
            "aspartate_aminotransferase", "gamma_glutamyltransferase", "creatinine", "glucose",
            "total_protein","albumin", "testosterone"))

# Merge datasets ----
cohort <- mergeAllCovariates(read_csv(paste0(dir_data, "Results/",cohort_name,"_cohort.csv")),
                             baselineCharacteristics,
                             biomarkers) |>
  select(-c("eid", "asthma","liver_disease_mild"))

data_biomarkers <- cohort |> select(any_of(colnames(biomarkers)), state)

# Splines ----
resultsNonLinear <- tibble(
  "risk_factor" = as.character(),
  "prob" = as.double(),
  "lower" = as.double(),
  "upper" = as.double(),
  "pval" = as.double(),
  "significance" = as.double(),
  "name" = as.character()
)

for(name in colnames(data_biomarkers)[colnames(data_biomarkers) != "state"]){
  n <- name

  x <- data_biomarkers |>
    rename("risk_factor" = n) |>
    select("state", "risk_factor") |>
    arrange(risk_factor)

  fitOne <- glm(state ~ risk_factor, family = binomial(link="logit"), data=x)
  fitTwo <- glm(state ~ ns(risk_factor, 4), family = binomial(link="logit"), data = x)

  predictions <- predict(fitTwo, newdata =  x, se.fit=TRUE, type = "link", level = 0.95)

  resultsNonLinear <- resultsNonLinear |>
    rbind(
      tibble(risk_factor = x$risk_factor,
             prob = plogis(predictions$fit),
             lower = plogis(predictions$fit - 1.96 * predictions$se.fit),
             upper = plogis(predictions$fit + 1.96 * predictions$se.fit),
             pval = anova(fitOne, fitTwo, test = "Chisq")$`Pr(>Chi)`[2],
             significance = if_else(anova(fitOne, fitTwo, test = "Chisq")$`Pr(>Chi)`[2] < 0.05, 1, 0),
             name = n)
    )
}

if(cohort_name == "6-Pacs"){tit <- "B) PACS linearity assessment"}else{tit <- "A) Long COVID linearity assessment"}
p <- resultsNonLinear |>
  left_join(
    tibbleNameBiomarkers() |>
      full_join(tibbleNameComorbidities(), by = c("risk_factor", "risk_factor1")) |>
      rename("name" = "risk_factor"),
    by = "name"
  ) |>
  mutate(significance = if_else(significance == 1, "< 0.05", ">= 0.05")) |>
  ggplot(aes(x = risk_factor, y = prob, color = significance, fill = significance)) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.2) +
  xlab(str_to_sentence(n)) +
  xlab(" ") +
  ylab("Probability of having the outcome") +
  theme_bw(base_family = "Calibri") +
  facet_wrap(vars(risk_factor1), scales = "free") +
  labs(fill = "P-Value test for non-linear model:",
       color = "P-Value test for non-linear model:",
       title = tit)

ggsave(plot = p, filename = paste0(dir_results, "9-LinearityAnalysis/",cohort_name,"_plot.png"), width = 35, height = 20, unit = "cm", dpi = 600)

write.csv(resultsNonLinear, file = paste0(dir_results,"9-LinearityAnalysis/", cohort_name, "_list.csv"), row.names = FALSE)

rm(list = c("baselineCharacteristics", "biomarkers", "cohort", "data_biomarkers", "fitOne", "fitTwo", "p", "predictions", "resultsNonLinear","x",
            "n", "name"))

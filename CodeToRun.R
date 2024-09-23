rm(list = ls())
library(here)
library(ggplot2)
library(tidyr)
library(lubridate)
library(forcats)
library(DiagrammeR)
library(glmnet)
library(mice)
library(readxl)
library(readr)
library(stringr)
library(flextable)
library(fastDummies)
library(glue)
library(officer)
library(ftExtra)
library(dplyr)
library(splines)
library(forestploter)
library(grid)
library(car)

source(here("Study","Functions.R"))
dir_data    <- "D:/Projects/LongCovidAndPasc/"
dir_results <- "D:/Projects/LongCovidAndPasc/Results/"

# Prepare data
source(here("Study","1-LoadCovariates.R"))

source(here("Study","2-ExploratoryAnalysisBeforeProcessing.R"))

source(here("Study","3-DataCleaning.R"))

source(here("Study","4-ExploratoryAnalysisAfterProcessing.R"))

# Long covid cohorts -------------
n_symptoms <- 1
washout_period <- 30
source(here("Study","5-CreateLongCovidCohorts.R"))

write.csv(longCovid_cohort, paste0(dir_results,"5-LongCovid_cohort.csv"), row.names = FALSE)
write.csv(attr(longCovid_cases_cohort,    "cohort_attrition"), paste0(dir_results,"5-Attrition_longCovid_cases.csv"))
write.csv(attr(longCovid_controls_cohort, "cohort_attrition"), paste0(dir_results,"5-Attrition_longCovid_controls.csv"))

rm(list = c("x", "x1", "x2", "x3", "longCovid_cases_cohort", "longCovid_controls_cohort",
            "longCovid_cohort", "covid19_result", "health_questionnaire"))

# PACS cohort -------------
n_symptoms <- 1
washout_period <- 30
source(here("Study","6-CreatePACSCohorts.R"))

write.csv(attr(pacs_controls_cohort, "cohort_attrition"), paste0(dir_results,"6-Attrition_pacs_controls.csv"))
write.csv(attr(pacs_cases_cohort, "cohort_attrition"),    paste0(dir_results,"6-Attrition_pacs_cases.csv"))
write.csv(pacs_cohort, paste0(dir_results,"6-Pacs_cohort.csv"), row.names = FALSE)

rm(list = c("t", "t1", "t2", "t3", "pacs_cases_cohort", "pacs_controls_cohort",
            "pacs_cohort", "covid19_result", "hes_data", "sequela_table"))

# Baseline characteristics of the cohorts ----
source(here("Study","7-BaselineCharacteristicsOfTheCohorts.R"))

# Correlation analysis ----
source(here("Study","8-CorrelationAnalysis.R"))

# Linearity analysis ----
cohort_name <- "5-LongCovid"
source(here("Study", "9-LinearityAnalysis.R"))
cohort_name <- "6-Pacs"
source(here("Study", "9-LinearityAnalysis.R"))

# Lasso regression ----
cohort_name <- "5-LongCovid"
source(here("Study","10-LassoRegression.R"))
cohort_name <- "6-Pacs"
source(here("Study","10-LassoRegression.R"))

# Confidence intervals crude analysis ----
analysis <- "main"
cohort_name <- "5-LongCovid"
source(here("Study","11-Linearity_Crude.R"))
cohort_name <- "6-Pacs"
source(here("Study","11-Linearity_Crude.R"))

# Confidence intervals adjusting ----
cohort_name <- "5-LongCovid"
source(here("Study","12-Linearity_Adjust.R"))
cohort_name <- "6-Pacs"
source(here("Study","12-Linearity_Adjust.R"))

# Plot main analysis ----
cohort_name <- "5-LongCovid"
source(here("Study","13-PlotLinearity.R"))
cohort_name <- "6-Pacs"
source(here("Study","13-PlotLinearity.R"))

# Secondary analysis - Removing some biomarkers ----
cohort_name <- "5-LongCovid"
source(here("Study","14-RemovingBiomarkers.R"))
cohort_name <- "6-Pacs"
source(here("Study","14-RemovingBiomarkers.R"))

# Third analysis - Sex stratified ----
cohort_name <- "5-LongCovid"
source(here("Study","15-SexStratified.R"))
cohort_name <- "6-Pacs"
source(here("Study","15-SexStratified.R"))

# Forth analysis - Bonferroni correction ----
library(api2lm)
cohort_name <- "5-LongCovid"
source(here("Study","16-Linearity_Adjust_Bonferroni.R"))
cohort_name <- "6-Pacs"
source(here("Study","16-Linearity_Adjust_Bonferroni.R"))

# Secondary analysis -----
n_symptoms <- 3
washout_period <- 30
source(here("Study","5-CreateLongCovidCohorts.R"))
write.csv(longCovid_cohort, paste0(dir_results,"SecondaryAnalysis/5-LongCovid_cohort.csv"), row.names = FALSE)
write.csv(attr(longCovid_cases_cohort,    "cohort_attrition"), paste0(dir_results,"SecondaryAnalysis/5-Attrition_longCovid_cases.csv"))
write.csv(attr(longCovid_controls_cohort, "cohort_attrition"), paste0(dir_results,"SecondaryAnalysis/5-Attrition_longCovid_controls.csv"))

rm(list = c("x", "x1", "x2", "x3", "longCovid_cases_cohort", "longCovid_controls_cohort",
            "longCovid_cohort", "covid19_result", "health_questionnaire"))

analysis <- "secondary"
cohort_name <- "5-LongCovid"
source(here("Study","11-Linearity_Crude.R"))

cohort_name <- "5-LongCovid"
source(here("Study","12-Linearity_Adjust.R"))

cohort_name <- "5-LongCovid"
source(here("Study","13-PlotLinearity.R"))


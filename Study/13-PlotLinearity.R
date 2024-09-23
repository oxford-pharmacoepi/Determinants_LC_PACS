windowsFonts("Calibri" = windowsFont("Calibri"))

get_scale <- function(plot,
                      width_wanted,
                      height_wanted,
                      unit = "in"){
  h <- convertHeight(sum(plot$heights), unit, TRUE)
  w <- convertWidth(sum(plot$widths), unit, TRUE)
  max(c(w/width_wanted,  h/height_wanted))
}

path <- if_else(analysis == "secondary",
                paste0(dir_results,"SecondaryAnalysis/",cohort_name,"_LogisticResults_Crude.csv"),
                paste0(dir_results,cohort_name,"_LogisticResults_Crude.csv"))
crude <- tibbleNameComorbidities() |>
  mutate(Type = "Comorbidities") |>
  full_join(tibbleNameBiomarkers() |>
              mutate(Type = "Biomarkers"),
            by = c("risk_factor", "risk_factor1", "Type")) |>
  full_join(
    as_tibble(read.csv(path)) |>
      mutate(Variable = if_else(Variable == "-","",Variable)) |>
      rename("risk_factor" = "Risk.factor", "PValue" = "Pr...z..") |>
      select(-"X"),
    by = c("risk_factor", "Type")
  ) |>
  filter((!is.na(Variable)) | (is.na(Variable) & risk_factor == risk_factor1)) |>
  mutate(Variable = gsub(" - ","",Variable)) |>
  mutate(Variable = if_else(is.na(Variable), "", Variable)) |>
  mutate(Type = if_else(is.na(Type), "Baseline characteristics", Type)) |>
  mutate(Analysis = "Crude")

path <- if_else(analysis == "secondary",
                paste0(dir_results,"SecondaryAnalysis/",cohort_name,"_LogisticResults_Adjusted.csv"),
                paste0(dir_results,cohort_name,"_LogisticResults_Adjusted.csv"))
adjusted <- tibbleNameComorbidities() |>
  mutate(Type = "Comorbidities") |>
  full_join(tibbleNameBiomarkers() |>
              mutate(Type = "Biomarkers"),
            by = c("risk_factor", "risk_factor1", "Type")) |>
  full_join(
    as_tibble(read.csv(path)) |>
      rename("risk_factor" = "Risk.factor", "PValue" = "Pr...z..") |>
      select(-"X"),
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


l1 <- if_else(rep(cohort_name == "5-LongCovid",3), c(2,2,5),c(5,3,5))
col <- c("#5AAEB2", "#376274", "#BC5856", "#8E2723", "#F3CC45","#F09336")
nam <- c("Baseline characteristics", "Biomarkers", "Comorbidities")
nam1 <- c("A) Sociodemographic factors", "B) Biomarkers", "C) Comorbidities")
plist <- list()

for(i in 1:3){

  t1 <- crude |>
    rbind(adjusted) |>
    filter(Type == nam[i]) |>
    mutate(risk_factor1 = if_else(is.na(risk_factor1), str_to_sentence(gsub("_"," ",risk_factor)), risk_factor1))

  t_crude <- t1 |>
    filter(Analysis == "Crude" | is.na(Analysis)) |>
    rename("OR_crude" = "OR", "CI_LOW_crude" = "CI_LOW", "CI_HIGH_crude" = "CI_HIGH",
           "PValue_crude" = "PValue")
  t_adj   <- t1 |>
    filter(Analysis == "Adjusted" | is.na(Analysis)) |>
    rename("OR_adj" = "OR", "CI_LOW_adj" = "CI_LOW", "CI_HIGH_adj" = "CI_HIGH",
           "PValue_adj" = "PValue")

  if(cohort_name == "6-Pacs"){
    t_crude <- t_crude |> filter(risk_factor != "Liver")
    t_adj <- t_adj |> filter(risk_factor != "Liver")
  }
  rm("t1")

  dt <- t_crude |>
    select("risk_factor" = "risk_factor1", "Type", "Variable", ends_with("_crude")) |>
    left_join(
      t_adj |>
        select("risk_factor" = "risk_factor1", "Type", "Variable", ends_with("_adj")),
      by = c("risk_factor", "Type", "Variable")
    ) |>
    mutate(n = rev(row_number())) |>
    mutate(n = min(n, na.rm = TRUE) + row_number()-1, .by = c("risk_factor", "Type")) |>
    arrange(desc(n)) |>
    mutate(
      "OR" = if_else(
        is.na(OR_adj),
        paste0(round(OR_crude,2)," (", round(CI_LOW_crude,2),", ", round(CI_HIGH_crude,2),")"),
        paste0(round(OR_crude,2)," (", round(CI_LOW_crude,2),", ", round(CI_HIGH_crude,2),")\n",
               round(OR_adj,2)," (", round(CI_LOW_adj,2),", ", round(CI_HIGH_adj,2),")")
      ),
      "Pval" = if_else(
        is.na(OR_adj),
        if_else(PValue_crude < 0.05, "<0.05",paste0(formatC(PValue_crude, format = "e", 1))),
        paste0( if_else(PValue_crude < 0.05, "<0.05",paste0(formatC(PValue_crude, format = "e", 1))),"\n",
                if_else(PValue_adj < 0.05, "<0.05",paste0(formatC(PValue_adj, format = "e", 1))))
      )) |>
    mutate("OR" = if_else(is.na(OR_crude), "", OR),
           "Pval" = if_else(is.na(Pval), "", Pval)) |>
    arrange(rev(n))

  if(i == 1){
    dt <- dt |> arrange(risk_factor, Variable)
  }

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

  rm("t_crude","t_adj")

  if(i == 2 && cohort_name == "6-Pacs"){
    dt1 <- dt[1:28,]
    dt2 <- dt[29:nrow(dt),] |>
      # add_row() |> add_row() |> add_row() |>
      # add_row() |> add_row() |> add_row() |>
      mutate(across(is.character, function(x){if_else(is.na(x), "", x)})) |>
      rename("  " = "OR (95% CI)", "   " = "B) Biomarkers", "    " = "P-Value")

    plist[[2]] <- getForestPlot(dt1, dir_results, cohort_name, nam[i], l1[i], col[2+2*(i-1)], col[1+2*(i-1)], extra = 1)
    plist[[4]] <- getForestPlot(dt2, dir_results, cohort_name, nam[i], l1[i], col[2+2*(i-1)], col[1+2*(i-1)], extra = 2)

  }else{
    plist[[i]] <- getForestPlot(dt, dir_results, cohort_name, nam[i], l1[i],col[2+2*(i-1)], col[1+2*(i-1)])
  }
}

max_width <- unit.pmax(plist[[1]]$widths, plist[[3]]$widths, plist[["2.1"]]$widths, plist[["2.2"]]$widths)

# Set the same widths for both plots
i <- 1
for(p in plist){
  p$widths <- max_width
  if(i == 2 | i == 4){
    h <- 16
    w <- 11
  }else{
    h <- 7
    w = 8
    }
  p_sc <- get_scale(plot = p, width_wanted = w, height_wanted = h, unit = "in")
  ggplot2::ggsave(filename = if_else(analysis == "secondary",
                                     paste0(dir_results,"SecondaryAnalysis/",cohort_name,"_forestplot_main_",i,".png"),
                                     paste0(dir_results,"Figures/Main/",cohort_name,"_forestplot_main_",i,".png")),
                  plot = p,
                  dpi = 500,
                  units = "in",
                  width = w, height = h,
                  scale = p_sc
  )
  i <- i + 1
}

rm(list = c("crude","adjusted", "dataAnalysis","t","t1","i","l1","nam","path"))


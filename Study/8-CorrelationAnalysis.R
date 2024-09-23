# Load data ----
baselineCharacteristics <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/baselineCharacteristics.csv"))) |> select(-c("X"))
biomarkers <- as_tibble(read.csv(paste0(dir_results,"3-CleanData/biomarkers.csv")))  |> select(-c("X"))
comorbidities <- as_tibble(read.csv(paste0(dir_results, "3-CleanData/comorbidities.csv"))) |> select(-c("X"))

# Merge all data ----
data_clean <- baselineCharacteristics |>
  left_join(biomarkers, by = "eid") |>
  left_join(comorbidities, by = "eid") |>
  select(-c("year_of_birth","eid"))

data_correlation <- data_clean |>
  select(any_of(colnames(baselineCharacteristics)),
         any_of(sort(colnames(biomarkers))),
         any_of(sort(colnames(comorbidities))))

table_cor <- tibble(var1 = colnames(data_correlation)) |>
  cross_join(
    tibble(
      var2 = colnames(data_correlation)
    )
  ) |>
  mutate(value = NA)

for(i in 1:nrow(table_cor)){
  x <- cor(
    data_correlation[[table_cor[[i,"var1"]]]],
    data_correlation[[table_cor[[i,"var2"]]]],
    use = "complete.obs")

  table_cor <- table_cor |>
    mutate(value = if_else(var1 == table_cor[[i,"var1"]] &
                             var2 == table_cor[[i,"var2"]],
                           x,
                           value))
}

sociodemographics <- c("sex", "ethnic_background", "body_mass_index","smoking_status",
        "index_of_multiple_deprivation")

write.csv(table_cor, paste0(dir_results,"7-CorrelationAnalysis.csv"), row.names = FALSE)

table_cor |>
  filter(var1 %in% sociodemographics & var2 %in% sociodemographics) |>
  arrange(var1, var2) |>
  mutate(var1 = str_to_sentence(gsub("_"," ",var1))) |>
  mutate(var2 = str_to_sentence(gsub("_"," ",var2))) |>
  ggplot(aes(x = fct_inorder(var1), y = fct_inorder(var2), fill = value), color = "black") +
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "red",limits=c(-1, 1)) +
  xlab("") + ylab("") +
  labs(fill = "r", title = "A)")
ggsave(filename = paste0(dir_results,"Figures/Correlation_A.png"), width = 12, height = 10, units = "cm")


table_cor |>
  filter((var1 %in% (tibbleNameBiomarkers() |> pull("risk_factor"))) &
           (var2 %in% (tibbleNameBiomarkers() |> pull("risk_factor")))) |>
  left_join(tibbleNameBiomarkers() |> rename(var1 = risk_factor), by = "var1") |>
  left_join(tibbleNameBiomarkers() |> rename(var2 = risk_factor, risk_factor2 = risk_factor1), by = "var2") |>
  ggplot(aes(x = risk_factor1, y = risk_factor2, fill = value), color = "black") +
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "red",limits=c(-1, 1)) +
  xlab("") + ylab("") +
  labs(fill = "r", title = "B)")
ggsave(filename = paste0(dir_results,"Figures/Correlation_B.png"), width = 25, height = 20, units = "cm")


table_cor |>
  filter(var1 %in% (tibbleNameComorbidities() |> pull("risk_factor")) &
           var2 %in% (tibbleNameComorbidities() |> pull("risk_factor"))) |>
  left_join(tibbleNameComorbidities() |> rename(var1 = risk_factor), by = "var1") |>
  left_join(tibbleNameComorbidities() |> rename(var2 = risk_factor, risk_factor2 = risk_factor1), by = "var2") |>
  ggplot(aes(x = risk_factor1, y = risk_factor2, fill = value), color = "black") +
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "red",limits=c(-1, 1)) +
  xlab("") + ylab("") +
  labs(fill = "r", title = "C)")
ggsave(filename = paste0(dir_results,"Figures/Correlation_C.png"), width = 25, height = 20, units = "cm")


rm(list = c("i","baselineCharacteristics","biomarkers","comorbidities",
            "data_clean","data_correlation","table_cor","x"))

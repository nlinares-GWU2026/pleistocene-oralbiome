# Load and label Fellows Yates 2021 sample metadata (Neanderthal + Homo sapiens
# across Palaeolithic, Mesolithic, historical, and modern-day samples)
library(dplyr)

meta <- read.delim("data/raw/fellowsyates2021_sample_metadata.tsv", sep = "\t")

# Keep only real samples (drop lab contamination controls), and only Homo hosts
# (drop Gorilla, Pan, Alouatta comparison primates)
homo_meta <- meta %>%
  filter(Sample_or_Control == "Sample",
         Host_Common %in% c("Homo (Neanderthal)", "Homo (Modern Human)"))

# Translate each sample's Age label into our project's group scheme
homo_meta <- homo_meta %>%
  mutate(project_group = case_when(
    Host_Common == "Homo (Neanderthal)" ~ "Neanderthal",
    Age %in% c("Palaeolithic", "Epipalaeolithic", "Mesolithic", "Later_Stone_Age") ~ "Pre-agricultural sapiens",
    Age == "Historical" ~ "Post-agricultural/historical",
    Age == "ModernDay" ~ "Modern/living reference",
    TRUE ~ "Unclassified"
  ))

cat("Sample counts by project group:\n")
print(table(homo_meta$project_group))

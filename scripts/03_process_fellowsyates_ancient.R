# Load and label Fellows Yates 2021 sample metadata (Neanderthal + Homo sapiens
# across Palaeolithic, Mesolithic, historical, and modern-day samples), then
# connect those samples to their real genus-level abundance data.
library(dplyr)

# --- Metadata: which samples exist, and which project group each belongs to ---
meta <- read.delim("data/raw/fellowsyates2021_sample_metadata.tsv", sep = "\t")

homo_meta <- meta %>%
  filter(Sample_or_Control == "Sample",
         Host_Common %in% c("Homo (Neanderthal)", "Homo (Modern Human)"))

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

# --- Abundance data: connect each sample to its real genus-level counts ---
abundance <- read.delim("data/raw/fellowsyates2021_genus_abundance_nt.tsv", sep = "\t", check.names = FALSE)
abundance_cols <- colnames(abundance)[-1]

find_column <- function(sample_id) {
  hit <- abundance_cols[startsWith(abundance_cols, sample_id)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

homo_meta$abundance_column <- sapply(homo_meta$X.SampleID, find_column)
cat("Matched:", sum(!is.na(homo_meta$abundance_column)),
    "/ Unmatched:", sum(is.na(homo_meta$abundance_column)), "\n")

genus_names <- abundance[[1]]
ancient_counts <- abundance[, homo_meta$abundance_column]
colnames(ancient_counts) <- homo_meta$X.SampleID
rownames(ancient_counts) <- genus_names

# Convert to relative abundance (proportion per sample)
ancient_relabund <- sweep(as.matrix(ancient_counts), 2, colSums(ancient_counts), FUN = "/")

cat("Relative abundance sanity check (should all be ~1):\n")
print(summary(colSums(ancient_relabund)))

# Save results so this never needs recomputing from scratch
saveRDS(ancient_relabund, "data/processed/fellowsyates2021_genus_relabund.rds")
saveRDS(homo_meta, "data/processed/fellowsyates2021_sample_groups.rds")

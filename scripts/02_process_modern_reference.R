# Pull HMP 16S V3-V5 data and isolate the modern living-person comparison group
# (Supragingival Plaque = the living tissue that dental calculus forms from)

library(HMP16SData)
library(SummarizedExperiment)
library(dplyr)

v35 <- V35()

# Keep only Supragingival Plaque samples (tissue-matched to ancient dental calculus)
supra <- v35[, colData(v35)$HMP_BODY_SUBSITE == "Supragingival Plaque"]

# Collapse OTU-level counts up to genus level
counts <- assay(supra, "16SrRNA")
genus <- rowData(supra)$GENUS

df <- as.data.frame(counts)
df$GENUS <- ifelse(is.na(genus), "Unclassified", genus)

genus_counts <- df %>%
  group_by(GENUS) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Convert raw counts to relative abundance (proportion per sample)
genus_labels <- genus_counts$GENUS
count_matrix <- as.matrix(genus_counts[, -1])
relabund_matrix <- sweep(count_matrix, 2, colSums(count_matrix), FUN = "/")
rownames(relabund_matrix) <- genus_labels

# Save the processed result so we don't have to redo this every session
saveRDS(relabund_matrix, "data/processed/hmp_supragingival_genus_relabund.rds")

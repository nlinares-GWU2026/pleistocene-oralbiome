# Combine HMP (modern reference) and Fellows Yates 2021 (Neanderthal/ancient/
# historical/modern) genus-level relative abundance data into one
# comparison-ready table.
#
# Design note: HMP and Fellows Yates used different sequencing methods (16S vs.
# shotgun), and a Bray-Curtis check below shows a real batch effect even
# between their two independent modern-human groups. Fellows Yates 2021 is
# therefore treated as the primary time-axis comparison (Neanderthal ->
# pre-agricultural -> historical -> modern-day, all measured the same way),
# while HMP serves as a secondary, independently-collected robustness check.

library(vegan)

hmp <- readRDS("data/processed/hmp_supragingival_genus_relabund.rds")
ancient <- readRDS("data/processed/fellowsyates2021_genus_relabund.rds")
homo_meta <- readRDS("data/processed/fellowsyates2021_sample_groups.rds")

hmp_genera <- rownames(hmp)
ancient_genera <- rownames(ancient)
shared_genera <- intersect(hmp_genera, ancient_genera)

cat("Shared genera used for comparison:", length(shared_genera), "\n")
cat("(HMP total:", length(hmp_genera), "| Fellows Yates total:", length(ancient_genera), ")\n")

hmp_shared <- hmp[shared_genera, ]
ancient_shared <- ancient[shared_genera, ]
combined <- cbind(hmp_shared, ancient_shared)

sample_info <- data.frame(
  sample_id = colnames(combined),
  project_group = c(rep("Modern/living reference (HMP)", ncol(hmp_shared)),
                     homo_meta$project_group[match(colnames(ancient_shared), homo_meta$X.SampleID)]),
  source_dataset = c(rep("HMP16SData", ncol(hmp_shared)),
                      rep("FellowsYates2021", ncol(ancient_shared))),
  analysis_role = c(rep("secondary_check", ncol(hmp_shared)),
                     rep("primary", ncol(ancient_shared)))
)

cat("\nSample counts by project group:\n")
print(table(sample_info$project_group))

# --- Batch effect diagnostic: do the two independent modern-human groups agree? ---
modern_only <- sample_info$project_group %in% c("Modern/living reference", "Modern/living reference (HMP)")
modern_dist_matrix <- as.matrix(vegdist(t(combined[, modern_only]), method = "bray"))

hmp_ids <- sample_info$sample_id[sample_info$project_group == "Modern/living reference (HMP)"]
fy_ids  <- sample_info$sample_id[sample_info$project_group == "Modern/living reference"]

within_hmp <- mean(modern_dist_matrix[hmp_ids, hmp_ids][lower.tri(modern_dist_matrix[hmp_ids, hmp_ids])])
between_groups <- mean(modern_dist_matrix[hmp_ids, fy_ids])

cat("\nBatch-effect check (Bray-Curtis distance):\n")
cat("  Within HMP modern samples:", round(within_hmp, 3), "\n")
cat("  Between HMP and Fellows-Yates modern samples:", round(between_groups, 3), "\n")

# Save the combined, comparison-ready data
saveRDS(combined, "data/processed/combined_genus_relabund.rds")
saveRDS(sample_info, "data/processed/combined_sample_info.rds")

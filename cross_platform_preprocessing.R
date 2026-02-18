#########################################
# Cross-Platform Integration: Microarray to RNA-seq
# Preprocessing TARGET RNA-seq data for integration with microarray data
#########################################

library(limma)
library(sva)  # for ComBat
library(edgeR)  # for RNA-seq normalization
library(dplyr)

# Create output directory for preprocessing results
preprocess_dir <- "Limma/preprocessing"
if (!dir.exists(preprocess_dir)) {
  dir.create(preprocess_dir, recursive = TRUE)
  cat("Created preprocessing directory:", preprocess_dir, "\n")
}

cat("\n=== CROSS-PLATFORM INTEGRATION PREPROCESSING ===\n")
cat("Microarray (Training) -> RNA-seq (TARGET Validation)\n\n")

#########################################
# STEP 1: Load and prepare TARGET RNA-seq data
#########################################

cat("STEP 1: Loading TARGET RNA-seq data...\n")
target_raw <- read.csv("Target_data.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Extract gene names and expression data
gene_names_target <- target_raw[[1]]
target_expression <- target_raw[, 2:(ncol(target_raw))]
rownames(target_expression) <- gene_names_target

cat("TARGET data dimensions:", dim(target_expression), "\n")
cat("Number of samples:", ncol(target_expression), "\n")
cat("Number of genes:", nrow(target_expression), "\n")

# Convert to numeric matrix
target_expression <- apply(target_expression, 2, function(x) as.numeric(as.character(x)))
rownames(target_expression) <- gene_names_target

# Remove genes with all zeros or NAs
keep_genes <- rowSums(target_expression, na.rm = TRUE) > 0
target_expression <- target_expression[keep_genes, ]
cat("After removing zero-expression genes:", nrow(target_expression), "genes\n")

#########################################
# STEP 2: RNA-seq specific normalization (TMM + log2)
#########################################

cat("\nSTEP 2: Applying RNA-seq normalization (TMM + log2)...\n")

# Create DGEList object for edgeR normalization
dge <- DGEList(counts = target_expression)

# Calculate normalization factors (TMM)
dge <- calcNormFactors(dge, method = "TMM")
cat("TMM normalization factors calculated\n")

# Get normalized counts (CPM - Counts Per Million)
target_cpm <- cpm(dge, log = FALSE)

# Apply log2 transformation with offset
target_log2 <- log2(target_cpm + 1)

cat("RNA-seq data normalized: log2(CPM + 1)\n")
cat("Expression range: [", round(min(target_log2), 2), ",", round(max(target_log2), 2), "]\n")

# Save intermediate result
write.csv(target_log2, 
          file.path(preprocess_dir, "target_log2cpm_normalized.csv"),
          row.names = TRUE)

#########################################
# STEP 3: Load microarray data for comparison
#########################################

cat("\nSTEP 3: Loading microarray training data...\n")
sig_genes_file <- "Limma/significant_genes_expression_WT_vs_NT.csv"

if (!file.exists(sig_genes_file)) {
  stop("Microarray significant genes file not found. Please run the main limma analysis first.")
}

microarray_data <- read.csv(sig_genes_file, row.names = 1, check.names = FALSE)
cat("Microarray data dimensions:", dim(microarray_data), "\n")

#########################################
# STEP 4: Find common genes between platforms
#########################################

cat("\nSTEP 4: Identifying common genes...\n")

microarray_genes <- rownames(microarray_data)
rnaseq_genes <- rownames(target_log2)

common_genes <- intersect(microarray_genes, rnaseq_genes)
cat("Microarray genes:", length(microarray_genes), "\n")
cat("RNA-seq genes:", length(rnaseq_genes), "\n")
cat("Common genes:", length(common_genes), "\n")

missing_in_rnaseq <- setdiff(microarray_genes, rnaseq_genes)
cat("Missing in RNA-seq:", length(missing_in_rnaseq), "\n")

if (length(missing_in_rnaseq) > 0) {
  cat("First 10 missing genes:", paste(head(missing_in_rnaseq, 10), collapse = ", "), "\n")
  writeLines(missing_in_rnaseq, 
             file.path(preprocess_dir, "genes_missing_in_TARGET.txt"))
}

# Subset to common genes
microarray_common <- microarray_data[common_genes, , drop = FALSE]
rnaseq_common <- target_log2[common_genes, , drop = FALSE]

cat("Data subset to", length(common_genes), "common genes\n")

#########################################
# STEP 5: Check distribution before batch correction
#########################################

cat("\nSTEP 5: Checking data distributions...\n")

# Calculate statistics
microarray_stats <- data.frame(
  Mean = mean(as.matrix(microarray_common), na.rm = TRUE),
  Median = median(as.matrix(microarray_common), na.rm = TRUE),
  SD = sd(as.matrix(microarray_common), na.rm = TRUE),
  Min = min(as.matrix(microarray_common), na.rm = TRUE),
  Max = max(as.matrix(microarray_common), na.rm = TRUE)
)

rnaseq_stats <- data.frame(
  Mean = mean(as.matrix(rnaseq_common), na.rm = TRUE),
  Median = median(as.matrix(rnaseq_common), na.rm = TRUE),
  SD = sd(as.matrix(rnaseq_common), na.rm = TRUE),
  Min = min(as.matrix(rnaseq_common), na.rm = TRUE),
  Max = max(as.matrix(rnaseq_common), na.rm = TRUE)
)

cat("\nMicroarray statistics:\n")
print(microarray_stats)
cat("\nRNA-seq statistics:\n")
print(rnaseq_stats)

# Create distribution plots
library(ggplot2)
library(reshape2)

# Prepare data for plotting
micro_melt <- melt(as.matrix(microarray_common))
micro_melt$Platform <- "Microarray"

rnaseq_melt <- melt(as.matrix(rnaseq_common))
rnaseq_melt$Platform <- "RNA-seq"

combined_melt <- rbind(micro_melt, rnaseq_melt)

# Plot distributions before correction
p1 <- ggplot(combined_melt, aes(x = value, fill = Platform)) +
  geom_density(alpha = 0.5) +
  labs(title = "Expression Distribution Before Batch Correction",
       x = "Expression Level", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Microarray" = "blue", "RNA-seq" = "red"))

ggsave(file.path(preprocess_dir, "distribution_before_correction.png"), 
       p1, width = 10, height = 6, dpi = 300)

#########################################
# STEP 6: Quantile normalization for cross-platform
#########################################

cat("\nSTEP 6: Applying quantile normalization for cross-platform integration...\n")

# Combine data
combined_data <- cbind(microarray_common, rnaseq_common)

# Create batch indicator
batch <- c(rep("Microarray", ncol(microarray_common)), 
           rep("RNAseq", ncol(rnaseq_common)))

cat("Total samples:", ncol(combined_data), "\n")
cat("Microarray samples:", sum(batch == "Microarray"), "\n")
cat("RNA-seq samples:", sum(batch == "RNAseq"), "\n")

# Apply quantile normalization across both platforms
combined_quantile <- normalizeBetweenArrays(combined_data, method = "quantile")

cat("Quantile normalization completed\n")

# Split back
microarray_quantile <- combined_quantile[, batch == "Microarray"]
rnaseq_quantile <- combined_quantile[, batch == "RNAseq"]

#########################################
# STEP 7: ComBat batch effect correction
#########################################

cat("\nSTEP 7: Applying ComBat batch effect correction...\n")

# Prepare batch variable
batch_factor <- factor(batch)

# Apply ComBat
# Note: We don't have covariates, so we use NULL
combat_data <- ComBat(dat = combined_quantile,
                      batch = batch_factor,
                      mod = NULL,
                      par.prior = TRUE,
                      prior.plots = FALSE)

cat("ComBat correction completed\n")

# Split back
microarray_combat <- combat_data[, batch == "Microarray"]
rnaseq_combat <- combat_data[, batch == "RNAseq"]

# Statistics after ComBat
microarray_combat_stats <- data.frame(
  Mean = mean(as.matrix(microarray_combat), na.rm = TRUE),
  Median = median(as.matrix(microarray_combat), na.rm = TRUE),
  SD = sd(as.matrix(microarray_combat), na.rm = TRUE),
  Min = min(as.matrix(microarray_combat), na.rm = TRUE),
  Max = max(as.matrix(microarray_combat), na.rm = TRUE)
)

rnaseq_combat_stats <- data.frame(
  Mean = mean(as.matrix(rnaseq_combat), na.rm = TRUE),
  Median = median(as.matrix(rnaseq_combat), na.rm = TRUE),
  SD = sd(as.matrix(rnaseq_combat), na.rm = TRUE),
  Min = min(as.matrix(rnaseq_combat), na.rm = TRUE),
  Max = max(as.matrix(rnaseq_combat), na.rm = TRUE)
)

cat("\nAfter ComBat - Microarray statistics:\n")
print(microarray_combat_stats)
cat("\nAfter ComBat - RNA-seq statistics:\n")
print(rnaseq_combat_stats)

#########################################
# STEP 8: Visualize after batch correction
#########################################

cat("\nSTEP 8: Creating visualization of batch correction...\n")

# Prepare data for plotting
micro_combat_melt <- melt(as.matrix(microarray_combat))
micro_combat_melt$Platform <- "Microarray"

rnaseq_combat_melt <- melt(as.matrix(rnaseq_combat))
rnaseq_combat_melt$Platform <- "RNA-seq"

combined_combat_melt <- rbind(micro_combat_melt, rnaseq_combat_melt)

# Plot distributions after correction
p2 <- ggplot(combined_combat_melt, aes(x = value, fill = Platform)) +
  geom_density(alpha = 0.5) +
  labs(title = "Expression Distribution After ComBat Correction",
       x = "Expression Level", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Microarray" = "blue", "RNA-seq" = "red"))

ggsave(file.path(preprocess_dir, "distribution_after_combat.png"), 
       p2, width = 10, height = 6, dpi = 300)

# Create boxplot comparison
sample_means_before <- data.frame(
  Mean = c(colMeans(microarray_common), colMeans(rnaseq_common)),
  Platform = batch,
  Stage = "Before"
)

sample_means_after <- data.frame(
  Mean = c(colMeans(microarray_combat), colMeans(rnaseq_combat)),
  Platform = batch,
  Stage = "After ComBat"
)

sample_means_all <- rbind(sample_means_before, sample_means_after)

p3 <- ggplot(sample_means_all, aes(x = Platform, y = Mean, fill = Platform)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~Stage) +
  labs(title = "Sample Mean Expression: Before and After ComBat",
       y = "Mean Expression") +
  theme_minimal() +
  scale_fill_manual(values = c("Microarray" = "blue", "RNAseq" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(preprocess_dir, "sample_means_comparison.png"), 
       p3, width = 10, height = 6, dpi = 300)

#########################################
# STEP 9: PCA to check batch effects
#########################################

cat("\nSTEP 9: Running PCA analysis...\n")

# Before correction
pca_before <- prcomp(t(combined_data), scale. = TRUE)
pca_before_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Platform = batch
)

p4 <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA Before Batch Correction",
       x = paste0("PC1 (", round(summary(pca_before)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_before)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Microarray" = "blue", "RNAseq" = "red"))

ggsave(file.path(preprocess_dir, "pca_before_correction.png"), 
       p4, width = 8, height = 6, dpi = 300)

# After correction
pca_after <- prcomp(t(combat_data), scale. = TRUE)
pca_after_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Platform = batch
)

p5 <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Platform)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA After ComBat Correction",
       x = paste0("PC1 (", round(summary(pca_after)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_after)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  scale_color_manual(values = c("Microarray" = "blue", "RNAseq" = "red"))

ggsave(file.path(preprocess_dir, "pca_after_correction.png"), 
       p5, width = 8, height = 6, dpi = 300)

#########################################
# STEP 10: Save preprocessed TARGET data
#########################################

cat("\nSTEP 10: Saving preprocessed TARGET data...\n")

# Save the ComBat-corrected RNA-seq data (ready for validation)
write.csv(rnaseq_combat, 
          file.path("Limma", "TARGET_preprocessed_combat_corrected.csv"),
          row.names = TRUE)

# Save the full combat-corrected data (both platforms)
write.csv(combat_data,
          file.path(preprocess_dir, "combined_data_combat_corrected.csv"),
          row.names = TRUE)

# Save quantile-normalized version (alternative)
write.csv(rnaseq_quantile,
          file.path(preprocess_dir, "TARGET_quantile_normalized.csv"),
          row.names = TRUE)

# Create summary report
summary_report <- c(
  "=== CROSS-PLATFORM INTEGRATION SUMMARY ===",
  "",
  paste("Date:", Sys.time()),
  "",
  "INPUT DATA:",
  paste("- Microarray genes:", length(microarray_genes)),
  paste("- RNA-seq genes:", length(rnaseq_genes)),
  paste("- Common genes:", length(common_genes)),
  paste("- Missing in TARGET:", length(missing_in_rnaseq)),
  "",
  "NORMALIZATION STEPS:",
  "1. RNA-seq: TMM normalization + log2(CPM+1)",
  "2. Quantile normalization across platforms",
  "3. ComBat batch effect correction",
  "",
  "STATISTICS BEFORE CORRECTION:",
  paste("  Microarray - Mean:", round(microarray_stats$Mean, 3), 
        "SD:", round(microarray_stats$SD, 3)),
  paste("  RNA-seq - Mean:", round(rnaseq_stats$Mean, 3), 
        "SD:", round(rnaseq_stats$SD, 3)),
  "",
  "STATISTICS AFTER COMBAT:",
  paste("  Microarray - Mean:", round(microarray_combat_stats$Mean, 3), 
        "SD:", round(microarray_combat_stats$SD, 3)),
  paste("  RNA-seq - Mean:", round(rnaseq_combat_stats$Mean, 3), 
        "SD:", round(rnaseq_combat_stats$SD, 3)),
  "",
  "OUTPUT FILES:",
  paste("- TARGET_preprocessed_combat_corrected.csv:", nrow(rnaseq_combat), "genes,", 
        ncol(rnaseq_combat), "samples"),
  "- Distribution plots (before/after)",
  "- PCA plots (before/after)",
  "- Sample mean comparison plots",
  "",
  "RECOMMENDATION:",
  "Use TARGET_preprocessed_combat_corrected.csv for validation analysis.",
  "This data has been:",
  "  1. Normalized using RNA-seq specific methods",
  "  2. Harmonized with microarray data",
  "  3. Batch-corrected using ComBat",
  "",
  "NEXT STEPS:",
  "1. Run survival analysis on preprocessed TARGET data",
  "2. Validate gene signature on ComBat-corrected data",
  "3. Compare results with non-corrected data as sensitivity analysis"
)

writeLines(summary_report, 
           file.path(preprocess_dir, "preprocessing_summary.txt"))

cat("\n=== PREPROCESSING COMPLETE ===\n")
cat("Preprocessed TARGET data saved to: Limma/TARGET_preprocessed_combat_corrected.csv\n")
cat("Summary report saved to:", file.path(preprocess_dir, "preprocessing_summary.txt"), "\n")
cat("\nYou can now use the preprocessed data for validation analysis.\n")

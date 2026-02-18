#########################################
# DEG Analysis with limma - For Microarray Data
# WT vs NT comparison
#########################################

# 1. Install and load necessary packages
# If you haven't installed them, uncomment the following lines:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")
# install.packages("dplyr")
# install.packages("tibble")
# install.packages("ggplot2")
# install.packages("pheatmap")
# install.packages("reshape2")
# install.packages("patchwork")

library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(patchwork)

# Create output directory
output_dir <- "Limma"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
}

# --- User-provided file paths ---
expression_file <- "merged_expression_dataset.csv"
survival_file <- "merged_labels.csv"

# 2. Load the expression data
cat("Loading expression data from:", expression_file, "\n")
expr_raw <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)

# The first column is empty, second column has gene symbols, rest are samples
gene_symbols <- expr_raw[[2]]
combined_expression <- expr_raw[, 3:ncol(expr_raw)]
rownames(combined_expression) <- gene_symbols

# Remove the gene symbol column name prefix from sample names
# Change "GSE2712.GSE2712_GSM52463" to "GSE2712_GSM52463"
original_colnames <- colnames(combined_expression)
colnames(combined_expression) <- gsub(".*\\.", "", colnames(combined_expression))

# Clean up any quotes or whitespace
colnames(combined_expression) <- trimws(colnames(combined_expression))
colnames(combined_expression) <- gsub("^\"|\"$", "", colnames(combined_expression))

cat("\nSample name transformation:\n")
cat("Original (first 5):", head(original_colnames, 5), "\n")
cat("After transformation (first 5):", head(colnames(combined_expression), 5), "\n")

# Convert to numeric matrix
# Use apply to ensure proper numeric conversion
combined_expression_numeric <- apply(combined_expression, 2, function(x) as.numeric(as.character(x)))
rownames(combined_expression_numeric) <- rownames(combined_expression)
colnames(combined_expression_numeric) <- colnames(combined_expression)
combined_expression <- combined_expression_numeric

cat("Dimensions of combined expression dataset:", dim(combined_expression), "\n")
cat("First few sample names:", head(colnames(combined_expression)), "\n")
cat("Data type check:", class(combined_expression), "Mode:", mode(combined_expression), "\n")

# Check for non-numeric values
if (any(is.na(combined_expression))) {
  na_count <- sum(is.na(combined_expression))
  cat("Warning:", na_count, "NA values detected after numeric conversion\n")
}

# Handle NA values by setting to 0
combined_expression[is.na(combined_expression)] <- 0

# Verify it's numeric
cat("First gene, first 5 samples:\n")
print(combined_expression[1, 1:min(5, ncol(combined_expression))])
cat("Is numeric?:", is.numeric(combined_expression), "\n")

# 3. Load the label data
cat("\nLoading label data from:", survival_file, "\n")
combined_labels <- read.csv(survival_file, stringsAsFactors = FALSE)

# Clean up any whitespace or special characters from sample names
combined_labels$sample <- trimws(combined_labels$sample)
combined_labels$sample <- gsub("^\"|\"$", "", combined_labels$sample)  # Remove quotes if present

cat("Dimensions of combined labels dataset:", dim(combined_labels), "\n")
cat("First few labels:\n")
print(head(combined_labels))

# 4. Prepare metadata for limma
sample_id_column_name <- "sample"

# Check if the sample_id_column_name exists in combined_labels
if (!sample_id_column_name %in% colnames(combined_labels)) {
  stop(paste0("Error: Column '", sample_id_column_name, "' not found in combined_labels.csv"))
}

# Get sample names from expression data
expression_sample_names <- colnames(combined_expression)

cat("\nDEBUGGING SAMPLE NAME MATCHING:\n")
cat("First 5 expression sample names:", head(expression_sample_names, 5), "\n")
cat("First 5 label sample names:", head(combined_labels[[sample_id_column_name]], 5), "\n")

# Check if samples exist in labels
samples_in_labels <- expression_sample_names %in% combined_labels[[sample_id_column_name]]
cat("Number of expression samples found in labels:", sum(samples_in_labels), "out of", length(expression_sample_names), "\n")

if (sum(samples_in_labels) == 0) {
  stop("ERROR: No matching samples found between expression data and labels file. Please check sample name formats.")
}

# Filter combined_labels to only include samples in expression data
coldata <- combined_labels %>%
  filter(!!sym(sample_id_column_name) %in% expression_sample_names)

cat("Filtered coldata has", nrow(coldata), "samples\n")

# Set row names
rownames(coldata) <- coldata[[sample_id_column_name]]
coldata <- coldata %>% select(-!!sym(sample_id_column_name))

# Reorder coldata rows to match expression matrix columns
# Only keep samples that exist in both
common_samples <- intersect(expression_sample_names, rownames(coldata))
cat("Common samples between expression and labels:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("ERROR: No common samples found after filtering.")
}

# Subset and reorder both matrices
combined_expression <- combined_expression[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

cat("Final dimensions - Expression:", dim(combined_expression), "Coldata:", dim(coldata), "\n")

# Check if 'label' column exists in coldata
if (!"label" %in% colnames(coldata)) {
  stop("Error: Column 'label' not found in the processed coldata.")
}

# Convert 'label' to a factor (WT vs NT)
coldata$label <- factor(coldata$label, levels = c("NT", "WT"))

# Verify that column names of expression matrix match row names of coldata
if (!all(colnames(combined_expression) == rownames(coldata))) {
  cat("ERROR: Mismatch detected!\n")
  cat("First 5 expression columns:", head(colnames(combined_expression), 5), "\n")
  cat("First 5 coldata rows:", head(rownames(coldata), 5), "\n")
  stop("Sample names in expression matrix and coldata do not match or are not in the same order.")
} else {
  cat("\nSample names in expression matrix and coldata match and are in order.\n")
}

cat("Label distribution:\n")
print(table(coldata$label))

# 5. Pre-filtering: Remove genes with very low expression or low variance
# For microarray, filter based on expression level and variance
gene_means <- rowMeans(combined_expression)
gene_vars <- apply(combined_expression, 1, var)

# Keep genes with reasonable expression (above 25th percentile of mean) and variance
keep_mean <- gene_means > quantile(gene_means, 0.25)
keep_var <- gene_vars > 0  # Remove genes with zero variance

keep <- keep_mean & keep_var
filtered_expression <- combined_expression[keep, ]

cat("\nRemoved", nrow(combined_expression) - nrow(filtered_expression), "genes with low expression or zero variance.\n")
cat("Remaining genes for analysis:", nrow(filtered_expression), "\n")

# 6. Normalize the data using quantile normalization (standard for microarray)
cat("\nPerforming quantile normalization...\n")
normalized_expression <- normalizeBetweenArrays(filtered_expression, method = "quantile")
cat("Normalization complete.\n")

# 7. Set up the design matrix
cat("\nSetting up design matrix...\n")
design <- model.matrix(~ 0 + label, data = coldata)
colnames(design) <- levels(coldata$label)
cat("Design matrix:\n")
print(head(design))

# 8. Fit the linear model
cat("\nFitting linear model...\n")
fit <- lmFit(normalized_expression, design)

# 9. Set up contrasts (WT vs NT)
cat("Setting up contrasts for WT vs NT comparison...\n")
contrast_matrix <- makeContrasts(
  WT_vs_NT = WT - NT,
  levels = design
)
print(contrast_matrix)

# 10. Fit contrasts
cat("\nFitting contrasts...\n")
fit2 <- contrasts.fit(fit, contrast_matrix)

# 11. Apply empirical Bayes smoothing
cat("Applying empirical Bayes smoothing...\n")
fit2 <- eBayes(fit2)

# 12. Extract results
cat("\nExtracting results for WT vs NT comparison...\n")
results <- topTable(fit2, coef = "WT_vs_NT", number = Inf, sort.by = "P")

# Rename columns to match DESeq2 output format for consistency
results_formatted <- results %>%
  dplyr::rename(
    log2FoldChange = logFC,
    pvalue = P.Value,
    padj = adj.P.Val
  )

# Order by adjusted p-value
results_ordered <- results_formatted[order(results_formatted$padj), ]

# 13. Summarize results
cat("\n=== SUMMARY OF LIMMA RESULTS (WT vs NT) ===\n")
cat("Total genes analyzed:", nrow(results_ordered), "\n")
cat("Genes with padj < 0.05:", sum(results_ordered$padj < 0.05, na.rm = TRUE), "\n")
cat("Genes with padj < 0.05 & |log2FC| > 1:", 
    sum(results_ordered$padj < 0.05 & abs(results_ordered$log2FoldChange) > 1, na.rm = TRUE), "\n")

# 14. Display top differentially expressed genes
cat("\nTop 20 differentially expressed genes (by adjusted p-value):\n")
print(head(results_ordered[, c("log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B")], 20))

# Save the results to a CSV file
write.csv(results_ordered, 
          file = file.path(output_dir, "DEG_results_WT_vs_NT.csv"),
          row.names = TRUE)
cat("\nFull DEG results saved to", file.path(output_dir, "DEG_results_WT_vs_NT.csv"), "\n")

# Get significant genes
significant_genes <- subset(results_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nNumber of significant genes (padj < 0.05, |log2FoldChange| > 1):", nrow(significant_genes), "\n")

###########################
# Visualization - WT vs NT
###########################

# Install EnhancedVolcano if needed
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

# Save expression data for significant genes
if (nrow(significant_genes) > 0) {
  significant_expression <- normalized_expression[rownames(significant_genes), ]
  write.csv(significant_expression, 
            file = file.path(output_dir, "significant_genes_expression_WT_vs_NT.csv"),
            row.names = TRUE)
  cat("\nExpression data for significant genes saved to", 
      file.path(output_dir, "significant_genes_expression_WT_vs_NT.csv"), "\n")
  cat("Number of significant genes:", nrow(significant_genes), "\n")
} else {
  cat("\nNo significant genes found based on the criteria (padj < 0.05, |log2FoldChange| > 1).\n")
  cat("Adjusting criteria to padj < 0.05 for visualization...\n")
  significant_genes <- subset(results_ordered, padj < 0.05)
  if (nrow(significant_genes) > 0) {
    significant_expression <- normalized_expression[rownames(significant_genes), ]
    write.csv(significant_expression, 
              file = file.path(output_dir, "significant_genes_expression_WT_vs_NT.csv"),
              row.names = TRUE)
    cat("Saved", nrow(significant_genes), "genes with padj < 0.05\n")
  }
}

# --- Volcano Plot for WT vs NT ---
cat("\nGenerating Volcano Plot...\n")

# Prepare data for volcano plot
volcano_data <- results_ordered
rownames(volcano_data) <- make.unique(rownames(volcano_data))

volcano_plot <- EnhancedVolcano(
  volcano_data,
  lab = rownames(volcano_data),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano Plot: WT vs NT Differential Gene Expression',
  subtitle = 'Adjusted p-value < 0.05, |Log2 Fold Change| > 1',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 1.5,
  labSize = 3.0,
  colAlpha = 0.8,
  legendLabels = c('Not significant', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
  col = c('grey', 'forestgreen', 'royalblue', 'red2'),
  boxedLabels = FALSE,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~italic(adj.P.Val)),
  caption = paste0('Total genes = ', nrow(results_ordered), '\n',
                   'Up-regulated in WT = ', sum(results_ordered$log2FoldChange > 1 & results_ordered$padj < 0.05, na.rm = TRUE), '\n',
                   'Down-regulated in WT = ', sum(results_ordered$log2FoldChange < -1 & results_ordered$padj < 0.05, na.rm = TRUE)),
  axisLabSize = 12,
  titleLabSize = 14,
  subtitleLabSize = 10,
  captionLabSize = 8,
  legendLabSize = 10,
  legendIconSize = 4.0
)

# Save the volcano plot
ggsave(
  filename = file.path(output_dir, "volcano_plot_WT_vs_NT.tiff"),
  plot = volcano_plot,
  device = "tiff",
  width = 7,
  height = 8,
  dpi = 300
)
cat("Volcano plot saved as", file.path(output_dir, "volcano_plot_WT_vs_NT.tiff"), "\n")

# --- Heatmap of Top 20 Differentially Expressed Genes ---
cat("\nGenerating Heatmap for Top 20 Differentially Expressed Genes...\n")

# Get the top 20 genes based on adjusted p-value
top20_genes <- head(results_ordered, 20)
top20_gene_ids <- rownames(top20_genes)

# Use normalized expression for visualization
expression_for_heatmap <- normalized_expression[top20_gene_ids, ]

# Scale the expression data for better visualization (z-score by row/gene)
scaled_expression_for_heatmap <- t(scale(t(expression_for_heatmap)))

# Prepare annotation for columns (samples)
annotation_col <- data.frame(Label = coldata$label)
rownames(annotation_col) <- rownames(coldata)

# Define colors for the annotation
ann_colors = list(
  Label = c(WT = "orange", NT = "purple")
)

# Generate the heatmap
pheatmap(
  scaled_expression_for_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Top 20 Differentially Expressed Genes (WT vs NT)",
  fontsize_row = 8,
  fontsize_col = 6,
  fontsize = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = file.path(output_dir, "heatmap_top20_genes_WT_vs_NT.tiff"),
  width = 8,
  height = 10,
  dpi = 300
)
cat("Heatmap for top 20 genes saved as", file.path(output_dir, "heatmap_top20_genes_WT_vs_NT.tiff"), "\n")

# --- Boxplots for Top 10 Differentially Expressed Genes ---
cat("\nGenerating Boxplots for Top 10 Differentially Expressed Genes...\n")

# Get the top 10 genes based on adjusted p-value
top10_genes <- head(results_ordered, 10)
top10_gene_ids <- rownames(top10_genes)

# Extract normalized expression data
expression_for_boxplot <- normalized_expression[top10_gene_ids, ]

# Convert to a data frame and add gene names
expression_df <- as.data.frame(expression_for_boxplot) %>%
  rownames_to_column(var = "Gene")

# Melt the data frame to long format
melted_expression <- melt(expression_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression_Level")

# Merge with coldata to get the 'label' information
merged_data_for_boxplot <- left_join(melted_expression,
                                     coldata %>% rownames_to_column(var = "Sample"),
                                     by = "Sample")

# Create a list to store individual plots
plot_list <- list()

# Loop through each of the top 10 genes to create a boxplot
for (gene in top10_gene_ids) {
  p <- ggplot(subset(merged_data_for_boxplot, Gene == gene),
              aes(x = label, y = Expression_Level, fill = label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
    scale_fill_manual(values = c("WT" = "#FF9933", "NT" = "#009966")) +
    labs(title = gene, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
    )
  plot_list[[gene]] <- p
}

# Combine all plots into a single figure
combined_boxplot_plot <- wrap_plots(plot_list, ncol = 5, nrow = 2) +
  plot_annotation(
    title = 'Top 10 Differentially Expressed Genes (WT vs. NT)',
    caption = 'Expression Level (Normalized)',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.0)
    )
  ) & theme(axis.title.y = element_text(size = 10, angle = 90))

# Save the combined boxplot plot
ggsave(
  filename = file.path(output_dir, "boxplot_top10_genes_WT_vs_NT.tiff"),
  plot = combined_boxplot_plot,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300
)
cat("Boxplots for top 10 genes saved as", file.path(output_dir, "boxplot_top10_genes_WT_vs_NT.tiff"), "\n")

# --- Additional Summary Statistics ---
cat("\n=== ADDITIONAL SUMMARY ===\n")
cat("Total samples:", ncol(normalized_expression), "\n")
cat("WT samples:", sum(coldata$label == "WT"), "\n")
cat("NT samples:", sum(coldata$label == "NT"), "\n")
cat("Total genes after filtering:", nrow(normalized_expression), "\n")

# Recount significant genes with strict criteria
significant_genes_strict <- subset(results_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
cat("Significant DEGs (padj < 0.05 & |log2FC| > 1):", nrow(significant_genes_strict), "\n")
if (nrow(significant_genes_strict) > 0) {
  cat("Up-regulated in WT:", sum(significant_genes_strict$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("Down-regulated in WT:", sum(significant_genes_strict$log2FoldChange < 0, na.rm = TRUE), "\n")
}

############# ML version ###############################################################
cat("\n=== PREPARING DATA FOR MACHINE LEARNING ===\n")

# Check if significant genes file exists
sig_genes_file <- file.path(output_dir, "significant_genes_expression_WT_vs_NT.csv")
if (file.exists(sig_genes_file)) {
  # Load the CSV file
  data <- read.csv(sig_genes_file, row.names = 1)
  label <- read.csv("merged_labels.csv")
  
  # Transpose the data
  transposed_data <- t(data)
  
  # Convert the transposed matrix to a data frame
  transposed_data <- as.data.frame(transposed_data)
  
  # Add row names as a column
  transposed_data$SampleID <- rownames(transposed_data)
  
  # Reorder columns to have SampleID as the first column
  transposed_data <- transposed_data[, c(ncol(transposed_data), 1:(ncol(transposed_data)-1))]
  
  # Save the transposed data to a new CSV file
  write.csv(transposed_data, 
            file.path(output_dir, "transposed_significant_genes_expression_WT_vs_NT_ML.csv"), 
            row.names = FALSE)
  
  cat("Transposed data has been saved to '", 
      file.path(output_dir, "transposed_significant_genes_expression_WT_vs_NT_ML.csv"), "'\n")
} else {
  cat("Skipping ML data preparation - no significant genes found.\n")
}



##########################################################################################
################################## Test Dataset #########################################
#################### Slice TARGET GENES ################################################
##################### CROSS-PLATFORM PREPROCESSING FOR TARGET DATA ######################
cat("\n=== CROSS-PLATFORM PREPROCESSING ===\n")

# Check if we need to preprocess TARGET data
target_preprocessed_file <- file.path(output_dir, "TARGET_preprocessed_combat_corrected.csv")

if (file.exists("Target_data.csv") && file.exists(sig_genes_file)) {
  if (!file.exists(target_preprocessed_file)) {
    cat("TARGET data found but not yet preprocessed.\n")
    cat("Running cross-platform integration preprocessing...\n")
    cat("This will normalize and batch-correct TARGET RNA-seq data with microarray data.\n\n")
    
    # Run the preprocessing script
    source("cross_platform_preprocessing.R")
    
  } else {
    cat("Preprocessed TARGET data already exists.\n")
    cat("Using:", target_preprocessed_file, "\n")
  }
} else {
  cat("Skipping preprocessing - required files not found\n")
  if (!file.exists("Target_data.csv")) {
    cat("  - Target_data.csv not found\n")
  }
  if (!file.exists(sig_genes_file)) {
    cat("  - Significant genes file not found\n")
  }
}

##################### Slice TARGET GENES (from preprocessed data) #######################
cat("\n=== SLICING TARGET GENES ===\n")

# Use preprocessed TARGET data if available, otherwise use raw
if (file.exists(target_preprocessed_file)) {
  cat("Using preprocessed (ComBat-corrected) TARGET data\n")
  target_data_to_use <- target_preprocessed_file
} else if (file.exists("Target_data.csv")) {
  cat("WARNING: Using raw TARGET data without preprocessing\n")
  cat("For cross-platform analysis, preprocessing is strongly recommended!\n")
  target_data_to_use <- "Target_data.csv"
} else {
  cat("No TARGET data available\n")
  target_data_to_use <- NULL
}

# Check if significant genes file exists and TARGET data exists
if (!is.null(target_data_to_use) && file.exists(sig_genes_file)) {
  library(dplyr)
  
  # Read the significant genes file
  sig_genes <- read.csv(sig_genes_file, stringsAsFactors = FALSE)
  
  # Get gene names from row names
  sig_gene_names <- rownames(read.csv(sig_genes_file, row.names = 1))
  
  # Read the TARGET data (preprocessed or raw)
  cat("Reading TARGET data from:", target_data_to_use, "\n")
  
  if (target_data_to_use == target_preprocessed_file) {
    # Preprocessed data already has row names as genes
    rna_seq_data_matrix <- read.csv(target_data_to_use, row.names = 1, check.names = FALSE)
    target_gene_names <- rownames(rna_seq_data_matrix)
    
    # Convert back to data frame format for consistency
    rna_seq_data <- data.frame(GeneSymbol = target_gene_names, rna_seq_data_matrix)
    
  } else {
    # Raw data needs processing
    rna_seq_data <- read.csv(target_data_to_use, stringsAsFactors = FALSE, check.names = FALSE)
    rna_seq_data$GeneSymbol <- rna_seq_data[[1]]
  }
  
  # Check if the gene symbol column exists
  if (!"GeneSymbol" %in% colnames(rna_seq_data)) {
    cat("Warning: Column 'GeneSymbol' not found in TARGET data. Skipping TARGET gene slicing.\n")
  } else {
    # Check for missing genes FIRST
    missing_genes <- setdiff(sig_gene_names, rna_seq_data$GeneSymbol)
    
    if (length(missing_genes) > 0) {
      cat("\n*** NOTE: The following", length(missing_genes), "genes were not found in TARGET data ***\n")
      cat(paste(head(missing_genes, 10), collapse = ", "), "\n")
      if (length(missing_genes) > 10) cat("... and", length(missing_genes) - 10, "more\n")
      
      # Save missing genes list
      missing_genes_file <- file.path(output_dir, "missing_genes_in_TARGET.txt")
      writeLines(c(
        paste("Missing genes in Target_data.csv:", length(missing_genes)),
        "",
        "These genes from your significant DEG list are not present in the TARGET dataset:",
        "",
        missing_genes,
        "",
        "POSSIBLE REASONS:",
        "1. Different gene annotations between platforms",
        "2. Low expression genes filtered out in TARGET",
        "3. Gene name updates/aliases",
        "4. Genes not covered by RNA-seq",
        "",
        "SOLUTION IMPLEMENTED:",
        "Proceeding with available genes only (recommended approach)",
        paste("Coverage:", round((1 - length(missing_genes)/length(sig_gene_names))*100, 2), "%"),
        "",
        paste("Date:", Sys.time())
      ), missing_genes_file)
      cat("Missing genes list saved to:", missing_genes_file, "\n\n")
      
      # Proceed with available genes
      cat("PROCEEDING with genes available in TARGET dataset...\n")
      sig_gene_names_available <- sig_gene_names[!(sig_gene_names %in% missing_genes)]
      cat("Using", length(sig_gene_names_available), "out of", length(sig_gene_names), "significant genes\n")
      cat("Coverage:", round((length(sig_gene_names_available)/length(sig_gene_names))*100, 2), "%\n")
      
    } else {
      cat("All significant genes were found in TARGET data.\n")
      sig_gene_names_available <- sig_gene_names
    }
    
    # Filter rows where GeneSymbol matches the AVAILABLE significant genes
    filtered_data <- rna_seq_data %>% 
      filter(GeneSymbol %in% sig_gene_names_available)
    
    # Reorder filtered_data to match the order of sig_gene_names_available
    filtered_data <- filtered_data[match(sig_gene_names_available, filtered_data$GeneSymbol), ]
    
    # Remove any NA rows
    filtered_data <- filtered_data[!is.na(filtered_data$GeneSymbol), ]
    
    # Report the number of genes successfully sliced
    cat("Number of genes successfully sliced:", nrow(filtered_data), "\n")
    
    # Save the filtered and ordered data
    output_filename <- if (target_data_to_use == target_preprocessed_file) {
      "sliced_TARGET_preprocessed_genes_ordered.csv"
    } else {
      "sliced_TARGET_raw_genes_ordered.csv"
    }
    
    write.csv(filtered_data, 
              file.path(output_dir, output_filename), 
              row.names = FALSE)
    
    cat("Filtered data saved to:", file.path(output_dir, output_filename), "\n")
    
    # Also save the list of genes actually used
    genes_used_file <- file.path(output_dir, "genes_used_in_TARGET_analysis.txt")
    writeLines(c(
      "=== GENES USED IN TARGET VALIDATION ANALYSIS ===",
      "",
      paste("Total significant genes from discovery:", length(sig_gene_names)),
      paste("Genes found in TARGET:", length(sig_gene_names_available)),
      paste("Missing genes:", length(missing_genes)),
      paste("Coverage:", round((length(sig_gene_names_available)/length(sig_gene_names))*100, 2), "%"),
      "",
      paste("Data source:", ifelse(target_data_to_use == target_preprocessed_file, 
                                   "Preprocessed (ComBat-corrected)", "Raw")),
      "",
      "Genes included in analysis:",
      "",
      sig_gene_names_available
    ), genes_used_file)
    cat("List of genes used saved to:", genes_used_file, "\n")
    
    # Save info about preprocessing status
    if (target_data_to_use == target_preprocessed_file) {
      cat("\n*** IMPORTANT ***\n")
      cat("Using ComBat batch-corrected data for cross-platform validation\n")
      cat("This accounts for systematic differences between microarray and RNA-seq\n")
    }
  }
} else {
  cat("Skipping TARGET gene slicing - required files not found\n")
}

########################### Prepare TARGET Clinical Data
cat("\n=== PREPARING TARGET CLINICAL DATA ===\n")

target_sliced_file <- file.path(output_dir, "sliced_Target_data_genes_ordered.csv")
if (file.exists(target_sliced_file) && file.exists("TARGET-WT.clinical.csv")) {
  library(dplyr)
  
  # Read the expression data
  expr_data <- read.csv(target_sliced_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Extract sample names: columns from 2 to second last
  sample_names <- colnames(expr_data)[2:(ncol(expr_data) - 1)]
  
  # Read the clinical data
  clinical_data <- read.csv("TARGET-WT.clinical.csv", stringsAsFactors = FALSE)
  
  # Standardize sample names in clinical_data to match expression data format
  clinical_data$sample <- gsub("-", ".", clinical_data$sample)
  clinical_data$sample <- sub("A$", "", clinical_data$sample)
  clinical_data$sample <- sub("\\.([0-9]{2})A\\.", ".\\1.", clinical_data$sample)
  clinical_data$sample <- sub("\\.([0-9]{2})A$", ".\\1", clinical_data$sample)
  
  # Filter clinical_data to only include samples that exist in sample_names
  filtered_clinical <- clinical_data %>%
    filter(sample %in% sample_names)
  
  # Check for missing samples
  missing_samples <- setdiff(sample_names, clinical_data$sample)
  if (length(missing_samples) > 0) {
    cat("The following", length(missing_samples), "samples from expression data were not found in clinical_data:\n")
    cat(paste(head(missing_samples, 10), collapse = ", "), "\n")
    if (length(missing_samples) > 10) cat("... and", length(missing_samples) - 10, "more\n")
  } else {
    cat("All samples from expression data were found in clinical_data.\n")
  }
  
  # Report the number of samples successfully sliced
  cat("Number of samples successfully sliced:", nrow(filtered_clinical), "\n")
  
  # Save the filtered clinical data
  write.csv(filtered_clinical, 
            file.path(output_dir, "sliced_clinical_data_uniform.csv"), 
            row.names = FALSE)
  
  cat("Filtered clinical data has been saved to '", 
      file.path(output_dir, "sliced_clinical_data_uniform.csv"), "'\n")
  
  #################################### Keep order the same in both Target data and labels
  
  # Reorder clinical_data to match the order of sample_names
  ordered_clinical <- filtered_clinical[match(sample_names, filtered_clinical$sample), ]
  
  # Remove any NA rows
  ordered_clinical <- ordered_clinical[!is.na(ordered_clinical$sample), ]
  
  # Report the number of samples
  cat("Number of samples in ordered clinical data:", nrow(ordered_clinical), "\n")
  
  # Save the ordered clinical data
  write.csv(ordered_clinical, 
            file.path(output_dir, "ordered_sliced_clinical_data_uniform.csv"), 
            row.names = FALSE)
  
  cat("Ordered clinical data has been saved to '", 
      file.path(output_dir, "ordered_sliced_clinical_data_uniform.csv"), "'\n")
  
} else {
  if (!file.exists(target_sliced_file)) {
    cat("Skipping TARGET clinical data preparation - sliced_Target_data_genes_ordered.csv not found.\n")
  }
  if (!file.exists("TARGET-WT.clinical.csv")) {
    cat("Skipping TARGET clinical data preparation - TARGET-WT.clinical.csv not found.\n")
  }
}

##################################### Order Full phenotype
cat("\n=== ORDERING FULL PHENOTYPE DATA ===\n")

if (file.exists(target_sliced_file) && file.exists("All_clinical_data.csv")) {
  library(dplyr)
  
  # Read the expression data to get sample order
  expr_data <- read.csv(target_sliced_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Extract sample names
  sample_names <- colnames(expr_data)[2:(ncol(expr_data) - 1)]
  
  # Read the clinical data
  clinical_data <- read.csv("All_clinical_data.csv", stringsAsFactors = FALSE)
  
  # Reorder clinical_data to match the order of sample_names
  ordered_clinical <- clinical_data[match(sample_names, clinical_data$sample), ]
  
  # Remove any NA rows
  ordered_clinical <- ordered_clinical[!is.na(ordered_clinical$sample), ]
  
  # Check for missing samples
  missing_samples <- setdiff(sample_names, clinical_data$sample)
  if (length(missing_samples) > 0) {
    cat("The following", length(missing_samples), "samples from expression data were not found in clinical_data:\n")
    cat(paste(head(missing_samples, 10), collapse = ", "), "\n")
    if (length(missing_samples) > 10) cat("... and", length(missing_samples) - 10, "more\n")
  } else {
    cat("All samples from expression data were found in clinical_data.\n")
  }
  
  # Report the number of samples
  cat("Number of samples in ordered clinical data:", nrow(ordered_clinical), "\n")
  
  # Save the ordered clinical data
  write.csv(ordered_clinical, 
            file.path(output_dir, "ordered_sliced_All_clinical_data_uniform.csv"), 
            row.names = FALSE)
  
  cat("Ordered clinical data has been saved to '", 
      file.path(output_dir, "ordered_sliced_All_clinical_data_uniform.csv"), "'\n")
} else {
  if (!file.exists(target_sliced_file)) {
    cat("Skipping full phenotype ordering - sliced_Target_data_genes_ordered.csv not found.\n")
  }
  if (!file.exists("All_clinical_data.csv")) {
    cat("Skipping full phenotype ordering - All_clinical_data.csv not found.\n")
  }
}

# ########################################################################################
# # Uniform the Normalization
# ########################################################################################
cat("\n=== UNIFORM NORMALIZATION FOR TRAINING AND TEST DATA ===\n")

training_file <- file.path(output_dir, "significant_genes_expression_WT_vs_NT.csv")
test_file <- file.path(output_dir, "sliced_Target_data_genes_ordered.csv")

if (file.exists(training_file) && file.exists(test_file)) {
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  
  # Step 1: Load the training data
  training_data <- read.csv(training_file, row.names = 1, check.names = FALSE)
  
  # Transpose to have samples as rows, genes as columns
  train_df <- as.data.frame(t(training_data))
  
  # Step 2: Load the external test data
  test_data <- read.csv(test_file, row.names = 1, check.names = FALSE)
  
  # Remove the last column (redundant GeneSymbol)
  test_data <- test_data[, -ncol(test_data)]
  
  # Transpose to have samples as rows, genes as columns
  test_df <- as.data.frame(t(test_data))
  
  # Step 3: Ensure common genes
  common_genes <- intersect(colnames(train_df), colnames(test_df))
  train_df <- train_df[, common_genes]
  test_df <- test_df[, common_genes]
  cat("Number of common genes:", length(common_genes), "\n")
  
  # Function to plot histogram
  plot_hist <- function(data, title, filename) {
    flat_data <- as.vector(as.matrix(data))
    df_plot <- data.frame(Value = flat_data)
    p <- ggplot(df_plot, aes(x = Value)) +
      geom_histogram(bins = 50, fill = "blue", color = "black") +
      theme_minimal() +
      labs(title = title, x = "Expression Value", y = "Frequency")
    ggsave(file.path(output_dir, filename), p, width = 8, height = 6)
    print(p)
  }
  
  # Function to print summary
  print_summary <- function(data, name) {
    flat <- as.vector(as.matrix(data))
    cat("\n", name, " Summary:\n")
    cat("Min:", min(flat, na.rm = TRUE), "Max:", max(flat, na.rm = TRUE), "\n")
    cat("Mean:", mean(flat, na.rm = TRUE), "Std:", sd(flat, na.rm = TRUE), "\n")
    cat("Median:", median(flat, na.rm = TRUE), "\n")
  }
  
  # Before normalization
  print_summary(train_df, "Training before")
  plot_hist(train_df, "Training Data Distribution (Before Normalization)", "train_before_hist.png")
  print_summary(test_df, "Test before")
  plot_hist(test_df, "Test Data Distribution (Before Normalization)", "test_before_hist.png")
  
  # Step 4: Apply log2(x + 1) transformation
  log_train <- log2(train_df + 1)
  log_test <- log2(test_df + 1)
  
  # After log
  print_summary(log_train, "Training after log")
  plot_hist(log_train, "Training Data Distribution (After Log)", "train_after_log_hist.png")
  print_summary(log_test, "Test after log")
  plot_hist(log_test, "Test Data Distribution (After Log)", "test_after_log_hist.png")
  
  # Step 5: Apply z-score normalization using training parameters (per gene)
  gene_means <- colMeans(log_train, na.rm = TRUE)
  gene_sds <- apply(log_train, 2, sd, na.rm = TRUE)
  
  # Avoid division by zero
  gene_sds[gene_sds == 0] <- 1
  
  train_z <- sweep(sweep(log_train, 2, gene_means, "-"), 2, gene_sds, "/")
  test_z <- sweep(sweep(log_test, 2, gene_means, "-"), 2, gene_sds, "/")
  
  # After z-score
  print_summary(train_z, "Training after z-score")
  plot_hist(train_z, "Training Data Distribution (After Z-Score)", "train_after_z_hist.png")
  print_summary(test_z, "Test after z-score")
  plot_hist(test_z, "Test Data Distribution (After Z-Score)", "test_after_z_hist.png")
  
  # Step 6: Save the normalized datasets
  write.csv(train_z, 
            file.path(output_dir, "normalized_training_data.csv"), 
            row.names = TRUE)
  write.csv(test_z, 
            file.path(output_dir, "normalized_test_data.csv"), 
            row.names = TRUE)
  
  # Confirmation
  cat("\nNormalization complete. Normalized files saved.\n")
  cat("You can now use 'normalized_training_data.csv' for training/testing\n")
  cat("and 'normalized_test_data.csv' for external validation.\n")
} else {
  if (!file.exists(training_file)) {
    cat("Skipping normalization - training file not found.\n")
  }
  if (!file.exists(test_file)) {
    cat("Skipping normalization - test file not found.\n")
  }
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All results have been saved to the '", output_dir, "' folder.\n")
cat("\nGenerated files:\n")
list.files(output_dir, full.names = FALSE)

































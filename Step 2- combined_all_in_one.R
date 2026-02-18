# Load required libraries
library(dplyr)
library(readr)
library(tibble)

# List of mutual genes dataset files
expression_files <- c(
  "GSE2712_mutual_genes.csv",
  "GSE11024_mutual_genes.csv", 
  "GSE66405_mutual_genes.csv",
  "GSE73209_mutual_genes.csv"
)

# List of label files
label_files <- c(
  "GSE2712_label.csv",
  "GSE11024_label.csv", 
  "GSE66405_label.csv",
  "GSE73209_label.csv"
)

# Dataset names
dataset_names <- c("GSE2712", "GSE11024", "GSE66405", "GSE73209")

# Step 1: Merge Gene Expression Datasets
cat("=== MERGING GENE EXPRESSION DATASETS ===\n")

# Load all expression datasets
expression_datasets <- list()

for (i in 1:length(expression_files)) {
  cat("Loading expression data:", expression_files[i], "\n")
  
  # Load dataset
  dataset <- read.csv(expression_files[i], row.names = 1, check.names = FALSE)
  
  # Verify that genes are in the same order (should be alphabetical)
  if (i > 1) {
    previous_genes <- rownames(expression_datasets[[i-1]])
    current_genes <- rownames(dataset)
    
    if (!identical(previous_genes, current_genes)) {
      cat("Warning: Gene order differs between datasets. Sorting alphabetically...\n")
      dataset <- dataset[order(rownames(dataset)), ]
    }
  }
  
  # Add dataset prefix to sample names to avoid conflicts
  colnames(dataset) <- paste0(dataset_names[i], "_", colnames(dataset))
  
  expression_datasets[[dataset_names[i]]] <- dataset
  
  cat("  Genes:", nrow(dataset), "| Samples:", ncol(dataset), "\n")
}

# Verify all datasets have the same genes in the same order
all_genes <- lapply(expression_datasets, rownames)
gene_check <- sapply(all_genes, function(x) identical(x, all_genes[[1]]))

if (all(gene_check)) {
  cat("✓ All datasets have the same genes in the same order.\n")
} else {
  cat("⚠ Gene order mismatch. Reordering datasets...\n")
  # Reorder all datasets to match the first one
  reference_genes <- rownames(expression_datasets[[1]])
  for (i in 2:length(expression_datasets)) {
    expression_datasets[[i]] <- expression_datasets[[i]][reference_genes, ]
  }
}

# Merge expression datasets by columns (samples)
cat("\nMerging expression datasets...\n")
merged_expression <- do.call(cbind, expression_datasets)

cat("Merged expression dataset dimensions:", dim(merged_expression), "\n")
cat("Genes:", nrow(merged_expression), "| Samples:", ncol(merged_expression), "\n")

# Save merged expression dataset
write.csv(merged_expression, "merged_expression_dataset.csv")
cat("Merged expression data saved to: merged_expression_dataset.csv\n")

# Step 2: Merge Label Files
cat("\n=== MERGING LABEL FILES ===\n")

# Load all label files
label_datasets <- list()

for (i in 1:length(label_files)) {
  cat("Loading labels:", label_files[i], "\n")
  
  # Load label file (assuming it has at least sample names and labels)
  labels <- read.csv(label_files[i], stringsAsFactors = FALSE)
  
  # Check structure of label file
  cat("  Columns in", label_files[i], ":", paste(colnames(labels), collapse = ", "), "\n")
  cat("  Number of samples:", nrow(labels), "\n")
  
  # Add dataset identifier
  labels$dataset <- dataset_names[i]
  
  # Add dataset prefix to sample names to match expression data
  if ("sample" %in% colnames(labels)) {
    labels$sample <- paste0(dataset_names[i], "_", labels$sample)
  } else if (ncol(labels) > 0) {
    # Assume first column contains sample names
    sample_col <- colnames(labels)[1]
    labels[[sample_col]] <- paste0(dataset_names[i], "_", labels[[sample_col]])
  }
  
  label_datasets[[dataset_names[i]]] <- labels
}

# Merge label files by rows
cat("\nMerging label files...\n")
merged_labels <- do.call(rbind, label_datasets)

# Reset row names
rownames(merged_labels) <- NULL

cat("Merged labels dimensions:", dim(merged_labels), "\n")
cat("Total samples in merged labels:", nrow(merged_labels), "\n")

# Save merged labels
write.csv(merged_labels, "merged_labels.csv", row.names = FALSE)
cat("Merged labels saved to: merged_labels.csv\n")

# Step 3: Verify Consistency
cat("\n=== VERIFYING CONSISTENCY ===\n")

# Check if sample names match between expression data and labels
expression_samples <- colnames(merged_expression)

# Extract sample names from labels (try different possible column names)
label_sample_col <- NULL
possible_sample_cols <- c("sample", "Sample", "ID", "Id", "sample_id", "Sample_ID")

for (col in possible_sample_cols) {
  if (col %in% colnames(merged_labels)) {
    label_sample_col <- col
    break
  }
}

if (is.null(label_sample_col)) {
  # Use first column as sample names
  label_sample_col <- colnames(merged_labels)[1]
  cat("Using column", label_sample_col, "as sample identifier\n")
}

label_samples <- merged_labels[[label_sample_col]]

# Check for exact match
if (identical(sort(expression_samples), sort(label_samples))) {
  cat("✓ Perfect match between expression samples and label samples.\n")
} else {
  cat("⚠ Sample name mismatch detected.\n")
  cat("Expression samples:", length(expression_samples), "\n")
  cat("Label samples:", length(label_samples), "\n")
  cat("Intersection:", length(intersect(expression_samples, label_samples)), "\n")
  
  # Show samples that are missing
  missing_in_labels <- setdiff(expression_samples, label_samples)
  missing_in_expression <- setdiff(label_samples, expression_samples)
  
  if (length(missing_in_labels) > 0) {
    cat("Samples in expression data but missing in labels:", length(missing_in_labels), "\n")
  }
  if (length(missing_in_expression) > 0) {
    cat("Samples in labels but missing in expression data:", length(missing_in_expression), "\n")
  }
  
  # Reorder labels to match expression data order
  merged_labels <- merged_labels[match(expression_samples, label_samples), ]
  cat("Labels reordered to match expression data sample order.\n")
}

# Step 4: Create Sample Information Table
cat("\n=== CREATING SAMPLE INFORMATION ===\n")

sample_info <- data.frame(
  sample_id = colnames(merged_expression),
  dataset = sapply(strsplit(colnames(merged_expression), "_"), function(x) x[1]),
  original_sample_id = sapply(strsplit(colnames(merged_expression), "_"), function(x) paste(x[-1], collapse = "_"))
)

# Add information from merged labels
if (exists("merged_labels")) {
  # Merge sample info with labels
  sample_info <- merge(sample_info, merged_labels, 
                       by.x = "sample_id", by.y = label_sample_col, 
                       all.x = TRUE)
}

write.csv(sample_info, "sample_information.csv", row.names = FALSE)
cat("Sample information saved to: sample_information.csv\n")

# Step 5: Final Summary
cat("\n=== FINAL SUMMARY ===\n")
cat("Merged Expression Dataset:\n")
cat("  - Genes:", nrow(merged_expression), "\n")
cat("  - Samples:", ncol(merged_expression), "\n")
cat("  - Datasets:", length(unique(sample_info$dataset)), "\n")

cat("\nMerged Labels:\n")
cat("  - Total samples:", nrow(merged_labels), "\n")
cat("  - Columns:", paste(colnames(merged_labels), collapse = ", "), "\n")

# Count samples per dataset
sample_counts <- table(sample_info$dataset)
cat("\nSamples per dataset:\n")
for (dataset in names(sample_counts)) {
  cat("  -", dataset, ":", sample_counts[dataset], "samples\n")
}

# Check for missing values in expression data
missing_values <- sum(is.na(merged_expression))
cat("\nMissing values in expression data:", missing_values, "\n")
if (missing_values > 0) {
  cat("Percentage of missing values:", round(missing_values / (nrow(merged_expression) * ncol(merged_expression)) * 100, 4), "%\n")
}

cat("\n=== FILES CREATED ===\n")
cat("1. merged_expression_dataset.csv - Combined gene expression data\n")
cat("2. merged_labels.csv - Combined sample labels/metadata\n")
cat("3. sample_information.csv - Detailed sample information\n")

cat("\n=== MERGING COMPLETE ===\n")
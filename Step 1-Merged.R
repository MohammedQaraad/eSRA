# Load required libraries
library(dplyr)
library(readr)

# List of dataset files
dataset_files <- c(
  "GSE2712_gene_expression_symbols.csv",
  "GSE11024_gene_expression_symbols.csv", 
  "GSE66405_gene_expression_symbols.csv",
  "GSE73209_gene_expression_symbols.csv"
)

# Names for the datasets (without extension)
dataset_names <- c("GSE2712", "GSE11024", "GSE66405", "GSE73209")

# Step 1: Load all datasets and extract gene symbols
gene_lists <- list()
datasets <- list()

cat("Loading datasets and extracting gene symbols...\n")

for (i in 1:length(dataset_files)) {
  cat("Loading:", dataset_files[i], "\n")
  
  # Load dataset
  dataset <- read.csv(dataset_files[i], row.names = 1, check.names = FALSE)
  
  # Store the dataset
  datasets[[dataset_names[i]]] <- dataset
  
  # Extract gene symbols (row names)
  gene_symbols <- rownames(dataset)
  gene_lists[[dataset_names[i]]] <- gene_symbols
  
  cat("  Number of genes in", dataset_names[i], ":", length(gene_symbols), "\n")
}

# Step 2: Find mutual genes across all datasets
mutual_genes <- Reduce(intersect, gene_lists)

cat("\n=== MUTUAL GENES ANALYSIS ===\n")
cat("Number of mutual genes across all datasets:", length(mutual_genes), "\n")
cat("Percentage of genes that are mutual:", 
    round(length(mutual_genes) / min(sapply(gene_lists, length)) * 100, 2), "%\n")

# Display first 20 mutual genes
cat("\nFirst 20 mutual genes:\n")
print(head(mutual_genes, 20))

# Step 3: Create new datasets with only mutual genes
mutual_datasets <- list()

cat("\nCreating mutual genes datasets...\n")

for (name in dataset_names) {
  # Get the original dataset
  original_data <- datasets[[name]]
  
  # Filter to keep only mutual genes
  mutual_data <- original_data[rownames(original_data) %in% mutual_genes, ]
  
  # Sort genes alphabetically for consistency
  mutual_data <- mutual_data[order(rownames(mutual_data)), ]
  
  # Store the mutual dataset
  mutual_datasets[[name]] <- mutual_data
  
  # Save the new dataset
  output_filename <- paste0(name, "_mutual_genes.csv")
  write.csv(mutual_data, output_filename)
  
  cat("Saved:", output_filename, "with", nrow(mutual_data), "genes\n")
}

# Step 4: Verify that all mutual datasets have the same genes
cat("\n=== VERIFICATION ===\n")
for (name in dataset_names) {
  mutual_genes_in_dataset <- rownames(mutual_datasets[[name]])
  cat("Genes in", name, "mutual dataset:", length(mutual_genes_in_dataset), "\n")
}

# Check if all datasets have exactly the same genes
all_same_genes <- all(sapply(mutual_datasets, function(x) {
  identical(rownames(x), rownames(mutual_datasets[[1]]))
}))

if (all_same_genes) {
  cat("✓ SUCCESS: All mutual datasets contain exactly the same genes in the same order.\n")
} else {
  cat("⚠ WARNING: Mutual datasets have different gene sets or order.\n")
}

# Step 5: Create a summary table of dataset sizes
summary_table <- data.frame(
  Dataset = dataset_names,
  Original_Genes = sapply(datasets, nrow),
  Mutual_Genes = sapply(mutual_datasets, nrow),
  Samples = sapply(datasets, ncol),
  Retention_Rate = round(sapply(mutual_datasets, nrow) / sapply(datasets, nrow) * 100, 2)
)

cat("\n=== DATASET SUMMARY ===\n")
print(summary_table)

# Step 6: Save the mutual genes list for reference
write.csv(data.frame(Gene_Symbol = sort(mutual_genes)), 
          "mutual_genes_list.csv", row.names = FALSE)
cat("\nMutual genes list saved to: mutual_genes_list.csv\n")

# Step 7: Optional - Create a combined overview of expression values
cat("\n=== EXPRESSION VALUE SUMMARY ===\n")
for (name in dataset_names) {
  dataset <- mutual_datasets[[name]]
  cat(name, "mutual dataset expression summary:\n")
  cat("  Min:", round(min(dataset, na.rm = TRUE), 2), "\n")
  cat("  Mean:", round(mean(as.matrix(dataset), na.rm = TRUE), 2), "\n")
  cat("  Max:", round(max(dataset, na.rm = TRUE), 2), "\n")
  cat("  Missing values:", sum(is.na(dataset)), "\n\n")
}

# Step 8: Create a Venn diagram (optional - requires additional package)
if (require(VennDiagram)) {
  venn_plot <- venn.diagram(
    x = gene_lists,
    category.names = dataset_names,
    filename = NULL,
    output = TRUE,
    height = 3000,
    width = 3000,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = c("lightblue", "lightgreen", "lightpink", "lightyellow"),
    alpha = 0.5,
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135, -135),
    cat.dist = c(0.055, 0.055, 0.055, 0.055),
    cat.fontfamily = "sans"
  )
  
  # Save the Venn diagram
  pdf("mutual_genes_venn_diagram.pdf", width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
  cat("Venn diagram saved to: mutual_genes_venn_diagram.pdf\n")
} else {
  cat("Install 'VennDiagram' package for visual representation of gene overlaps.\n")
}

cat("\n=== PROCESSING COMPLETE ===\n")
cat("Mutual genes datasets created successfully!\n")
cat("Total mutual genes found:", length(mutual_genes), "\n")
cat("New files created:\n")
for (name in dataset_names) {
  cat("  - ", name, "_mutual_genes.csv\n")
}
cat("  - mutual_genes_list.csv\n")
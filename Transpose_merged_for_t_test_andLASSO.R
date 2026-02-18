library(DESeq2)
library(dplyr)
library(tibble)

# --- User-provided file paths ---
expression_file <- "merged_expression_dataset.csv"

# Read the expression data
expression_data <- read.csv(expression_file, row.names = 1, check.names = FALSE)

# Transpose the data (swap rows and columns)
transposed_data <- t(expression_data) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID")

# Save the transposed data as CSV
output_file <- gsub("\\.csv$", "_transposed.csv", expression_file)
write.csv(transposed_data, file = output_file, row.names = FALSE)

cat("Transposition completed!\n")
cat("Original file:", expression_file, "\n")
cat("Transposed file saved as:", output_file, "\n")
# EndoGenomics
# Script 01: Inspect GWAS summary statistics

# Use a personal R library
dir.create("~/R/library", recursive = TRUE, showWarnings = FALSE)

.libPaths(c("~/R/library", .libPaths()))

# Use a personal R library
dir.create("~/R/library", recursive = TRUE, showWarnings = FALSE)

.libPaths(c("~/R/library", .libPaths()))# Install data.table if it is not already installed
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}

library(data.table)

# Change this to your actual downloaded filename
input_file <- "data/raw/endometriosis_gwas.tsv.gz"

# Read the data
gwas <- fread(input_file)

# Display dimensions and column names
cat("Number of variants:", nrow(gwas), "\n")
cat("Number of columns:", ncol(gwas), "\n")
print(names(gwas))

# Preview the data
print(head(gwas))

# Check column types
str(gwas)

# Check missing values
missing_counts <- colSums(is.na(gwas))
print(missing_counts)

# Summarise p-values
summary(gwas$p_value)

# Count variants below the conventional significance threshold
significant <- sum(
  gwas$p_value < 5e-8,
  na.rm = TRUE
)

cat("Variants with p < 5e-8:", significant, "\n")

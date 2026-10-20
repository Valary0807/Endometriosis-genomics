# EndoGenomics
# Script 02: Manhattan plot

# Set up a personal R package library
personal_lib <- path.expand("~/R/library")

dir.create(
  personal_lib,
  recursive = TRUE,
  showWarnings = FALSE
)

.libPaths(c(personal_lib, .libPaths()))

# Load required packages
library(data.table)
library(ggplot2)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

library(data.table)
library(ggplot2)

input_file <- "data/raw/endometriosis_gwas.tsv.gz"

gwas <- fread(input_file)

# Keep rows with valid chromosome, position and p-value
gwas <- gwas[
  !is.na(chromosome) &
  !is.na(base_pair_location) &
  !is.na(p_value) &
  p_value >= 0 &
  p_value <= 1
]

# Standardise chromosome labels
gwas[, chromosome := sub("^chr", "", as.character(chromosome),
                         ignore.case = TRUE)]

# Keep autosomes 1–22 for the initial plot
gwas[, chr_num := suppressWarnings(as.integer(chromosome))]
gwas <- gwas[chr_num %in% 1:22]

# Remove invalid genomic positions
gwas <- gwas[base_pair_location > 0]

# Avoid infinite values if a p-value is exactly zero
gwas[, neg_log10_p := -log10(pmax(p_value, 1e-300))]

# Calculate cumulative genomic positions
chr_info <- gwas[
  ,
  .(chr_length = max(base_pair_location)),
  by = chr_num
][order(chr_num)]

chr_info[, offset := shift(cumsum(chr_length), fill = 0)]

gwas <- merge(gwas, chr_info, by = "chr_num")
gwas[, cumulative_position := base_pair_location + offset]

# Centre x-axis labels on each chromosome
axis_info <- gwas[
  ,
  .(center = (min(cumulative_position) +
              max(cumulative_position)) / 2),
  by = chr_num
][order(chr_num)]

# Create the plot
p <- ggplot(
  gwas,
  aes(x = cumulative_position, y = neg_log10_p)
) +
  geom_point(aes(color = factor(chr_num %% 2)), size = 0.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_x_continuous(
    breaks = axis_info$center,
    labels = axis_info$chr_num
  ) +
  labs(
    title = "GWAS of Endometriosis",
    subtitle = "Dashed line: conventional genome-wide significance threshold",
    x = "Chromosome",
    y = expression(-log[10](p)),
    color = "Chromosome group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p)

# Save the figure
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "figures/endometriosis_manhattan.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)

cat("Plot saved to figures/endometriosis_manhattan.png\n")

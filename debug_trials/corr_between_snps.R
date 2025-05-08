library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)
library(jsonlite)
library(parallel)
library(hrbrthemes) 
library(dplyr)
library(data.table)
library(reshape2)
library(tidyr) 

snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
genetic_variants <- snp_data[, grepl("chr", colnames(snp_data))] 
genetic_variants <- genetic_variants[, apply(genetic_variants, 2, var) != 0]
exclRare <- TRUE
ctrl_corr <- TRUE
if (exclRare) {
  genetic_variants <- genetic_variants[, apply(genetic_variants, 2, sum) > dim(genetic_variants)[1]*2*0.05]
  genetic_variants <- genetic_variants[, apply(genetic_variants, 2, sum) < dim(genetic_variants)[1]*2*0.95]
}
if (ctrl_corr) {
  corr_mat <- cor(genetic_variants)
  Sigma.distance = as.dist(1 - abs(corr_mat))
  fit = hclust(Sigma.distance, method="single")
  corr_max <- 0.75
  clusters = cutree(fit, h=1-corr_max)
  representative_snps <- numeric(length(unique(clusters)))
  
  # Find representative SNP for each cluster
  for (k in unique(clusters)) {
    cluster_indices <- which(clusters == k)
    if (length(cluster_indices) == 1) {
      representative_snps[k] <- cluster_indices[1]
    } else {
      # Calculate sum of absolute correlations within the cluster
      cluster_corr_sum <- colSums(abs(corr_mat)[cluster_indices, cluster_indices])
      # Find the SNP with the highest sum of correlations
      representative_snps[k] <- cluster_indices[which.max(cluster_corr_sum)]
    }
  }
  genetic_variants <- genetic_variants[,representative_snps]
}

dim(genetic_variants)

# Compute correlation matrix
corr_mat <- cor(genetic_variants)

# Convert to long format
corr_df <- reshape2::melt(corr_mat)

# Filter to only show correlations with absolute value > 0.75
subset_corr_df <- subset(corr_df, abs(value) > 0.9)

# Plot heatmap
ggplot(subset_corr_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Heatmap of Genetic Variant Correlations (|r| > 0.75)",
       fill = "Correlation")


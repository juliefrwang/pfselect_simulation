setwd("/Users/juliew/Project/Knockoff/PFSelect_copy")
devtools::load_all()

library(PFSelectcopy)
library(glue)
library(jsonlite)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
standardize <- FALSE
excl_rare <- TRUE
ctrl_corr <- FALSE # if true features are dropped
FDR <- 0.2
model <- "Lasso" # "Lasso" or "Grplasso" 
seed <- 300

# excl_rare <- as.logical(args[1])
# standardize <- as.logical(args[2])
# FDR <- as.numeric(args[3])

setwd("~/Project/Knockoff/run_real_data_julie") #TOCHANGE

snp_filepath <- "../Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 
snps <- snp_data[, grepl("chr", colnames(snp_data))]
snps <- snps[, apply(snps, 2, var) != 0]

if (excl_rare) {
  snps <- snps[, apply(snps, 2, sum) > dim(snps)[1]*2*0.05]
  snps <- snps[, apply(snps, 2, sum) < dim(snps)[1]*2*0.95]
}
if (ctrl_corr) {
  corr_mat <- cor(snps)
  Sigma.distance = as.dist(1 - abs(corr_mat))
  fit = hclust(Sigma.distance, method="single")
  corr_max <- 0.9
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
  snps <- snps[,representative_snps]
}
pcs <- snp_data[,c("PC1", "PC2", "PC3", "PC4")]
# pcs <- as.data.frame(scale(pcs))
AD <- snp_data[,c("AD")]

set.seed(seed)
snps_knockoff <- generate_knockoff_data(snps)
save_path <- glue("run_standardize_{standardize}_exclRare_{excl_rare}_ctrlCorr_{ctrl_corr}_0.9_fdr_{FDR}_model_{model}_seed_{seed}")
dir.create(save_path)
write.csv(snps_knockoff, file = glue("{save_path}/final_snp_knockoff_data.csv"), row.names = FALSE)

results <- get_importance_matrices(
  genetic_variants = snps,
  genetic_variants_knockoff = snps_knockoff,
  additional_covariates = NULL,
  Z = pcs, 
  y = AD,
  n_folds = 5,
  standardize = standardize,
  FDR_rate = FDR,
  model = model
)


write.csv(as.matrix(results$coefs), file = glue("{save_path}/lasso_coefs.csv"), row.names = TRUE)
write.csv(results$scaled_selection_matrix, file = glue("{save_path}/scaled_selection_matrix.csv"), row.names = FALSE)
write.csv(results$selection_matrix, file = glue("{save_path}/selection_matrix.csv"), row.names = FALSE)
write.csv(results$W_statistic_matrix, file = glue("{save_path}/W_statistic_matrix.csv"), row.names = FALSE)
write.csv(results$q_values, file = glue("{save_path}/q_values.csv"), row.names = FALSE)

selected_snp_names <- colnames(snps)[apply(results$selection_matrix, 2, any)]
file_connection <- file(glue("{save_path}/selected_snps.txt"), "w") # Open file for writing
cat(glue("Selected {length(selected_snp_names)} features based on the knockoff filter:"), file = file_connection, sep = "\n")
cat(selected_snp_names, file = file_connection, sep = "\n")
close(file_connection) # Close the file connection



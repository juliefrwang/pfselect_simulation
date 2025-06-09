library(PFSelectcopy)
library(glue)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(reshape2)
library(glmnet)
source("~/Project/Knockoff/run_real_data_julie/scripts/fig3_heatmap_sd.R")

exclRare <- TRUE
ctrl_corr <- FALSE
snp_data_subgroup <- "ALL" # "EUR"

snp_filepath <- "~/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv" # rows: samples, cols: PCs, EUR, ..., SNPs, Outcome(AD)
snp_data <- read.csv(snp_filepath)
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 

SNPs <- snp_data[, grepl("chr", colnames(snp_data))]
SNPs <- SNPs[, apply(SNPs, 2, var) != 0]

if (snp_data_subgroup == "EUR") {
  snp_data <- snp_data[snp_data['ans_label'] == "EUR",]
  SNPs <- snp_data[, grepl("chr", colnames(snp_data))]
  SNPs <- SNPs[, apply(SNPs, 2, var) != 0]
}

pcs <- snp_data[, c("PC1", "PC2", "PC3", "PC4")]
eur <- snp_data$EUR
non_eur <- 1 - eur
afr <- snp_data$AFR
amr <- snp_data$AMR
y <- snp_data$AD

if (exclRare) {
  SNPs <- SNPs[, apply(SNPs, 2, sum) > dim(SNPs)[1]*2*0.05]
  SNPs <- SNPs[, apply(SNPs, 2, sum) < dim(SNPs)[1]*2*0.95]
}

if (ctrl_corr) {
  corr_mat <- cor(SNPs)
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
  SNPs <- SNPs[,representative_snps]
  }


set.seed(300)

################  1. fit linear model separately for each snp ################ 

# Initialize vector to store p-values
p_values <- numeric(ncol(SNPs))

# Perform marginal analysis for each SNP
for (j in 1:ncol(SNPs)) {
  model <- lm(y ~ SNPs[, j] + pcs[,1] + pcs[,2] + pcs[,3] + pcs[,4])
  p_values[j] <- summary(model)$coefficients[2, 4] # Extract p-value for SNP_j
}
p_values_df <- data.frame(SNP = colnames(SNPs), P_Value = p_values)
write.csv(p_values_df, glue("marginal_test_exclRare_{exclRare}_ctrl_corr_{ctrl_corr}_{snp_data_subgroup}.csv"), row.names = FALSE)

# Create a data frame for plotting
p_values_df <- data.frame(SNP_name = paste0(1:length(p_values), ": ", colnames(SNPs)), SNP = 1:length(p_values), P_Value = p_values, NegLogP = -log10(p_values))
highlight_df <- subset(p_values_df, NegLogP > -log10(5e-5))
ggplot(p_values_df, aes(x = SNP, y = NegLogP)) +
  geom_point() + 
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(5e-5), linetype = "dashed", color = "orange") +
  annotate("text", x = 0.9 * length(p_values), y = -log10(5e-5), vjust = -0.3, label = "y = 5e-5", color = "orange") + 
  geom_text(data = highlight_df, aes(label = SNP), vjust = 0.5, hjust = -0.8 ,size = 3) +
  # Create a legend with fullname
  geom_point(data = highlight_df, 
             aes(color = as.factor(SNP)), size = 2) +
  scale_color_manual(values = setNames(rep("blue", nrow(highlight_df)), highlight_df$SNP),
                     labels = setNames(highlight_df$SNP_name, highlight_df$SNP)) +
  labs(x = "SNPs", y = expression(paste("-log"[10], "(p-value)")), color = expression(paste("SNPs with p-val < -log"[10], "(5e-5)"))) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  ggtitle(glue("marginal test of candidate SNPs\nexclRare_{exclRare}_ctrl_corr_{ctrl_corr}_{snp_data_subgroup}"))
ggsave(glue("marginal_test/marginal_exclRare_{exclRare}_ctrl_corr_{ctrl_corr}_{snp_data_subgroup}.png"), last_plot(), width = 10, height = 7.4)

################  2. fit linear models on all snps together ################  
df <- data.frame(SNPs, pcs, y)
model <- lm(y ~ ., data = df)
p_values_2 <- summary(model)$coefficients[2:151, 4] # Extract p-value for SNP_j

p_values_df <- data.frame(SNP_name = paste0(1:length(p_values_2), ": ", colnames(SNPs)), SNP = 1:length(p_values_2), P_Value = p_values_2, NegLogP = -log10(p_values_2))
highlight_df <- subset(p_values_df, NegLogP > -log10(5e-3))
ggplot(p_values_df, aes(x = SNP, y = NegLogP)) +
  geom_point() + 
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(5e-3), linetype = "dashed", color = "orange") +
  annotate("text", x = 0.9 * length(p_values_2), y = -log10(5e-3), vjust = -0.3, label = "y = 5e-3", color = "orange") + 
  # Create a legend with fullname
  geom_point(data = highlight_df, 
             aes(color = as.factor(SNP)), size = 2) +
  scale_color_manual(values = setNames(rep("blue", nrow(highlight_df)), highlight_df$SNP),
                     labels = setNames(highlight_df$SNP_name, highlight_df$SNP)) +
  labs(x = "SNPs", y = expression(paste("-log"[10], "(p-value)")), color = expression(paste("SNPs with p-val < -log"[10], "(5e-3)"))) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  ggtitle(glue("multiple linear regression of candidate SNPs\nexclRare_{exclRare}_ctrl_corr_{ctrl_corr}_{snp_data_subgroup}"))
ggsave(glue("marginal_test/mlR_exclRare_{exclRare}_ctrl_corr_{ctrl_corr}_{snp_data_subgroup}.png"), last_plot(), height = 5.6, width = 7.3)

################ 3.LASSO test. ################  
df <- data.frame(SNPs, pcs, y)
ncol_df <- dim(df)[2]
lasso_model <- cv.glmnet(x=as.matrix(df[,1:ncol_df-1]), y=as.matrix(df[,ncol_df]), alpha = 1)
coefs <- coef(lasso_model, s = "lambda.min")
nonzero_features <- rownames(coefs)[which(coefs != 0)]
print(nonzero_features)

################  4.LASSO-Knockoff with no interacting terms ################  
snps_knockoff <- generate_knockoff_data(SNPs)
results <- get_importance_matrices(
  genetic_variants = SNPs,
  genetic_variants_knockoff = snps_knockoff,
  additional_covariates = pcs,
  Z = numeric(dim(SNPs)[1]), 
  y = y,
  n_folds = 5,
  standardize = TRUE,
  FDR_rate = 0.2
)
selected_snp_names <- colnames(SNPs)[apply(results$selection_matrix, 2, any)]
colnames(results$scaled_selection_matrix) <- colnames(SNPs)
all_W_heterogeneity_plot <- plot_heatmap(
  scaled_selection_matrix = results$scaled_selection_matrix,
  selected_snps = selected_snp_names,
  ancestry = non_eur,
  ancestry_name = "non-EUR",
  category = "W_all",
  save_path = "~/Project/Knockoff/run_real_data_julie/marginal_test"
)

################  4.LASSO-Knockoff with pcs interacting terms ################  
snps_knockoff <- generate_knockoff_data(SNPs)
results <- get_importance_matrices(
  genetic_variants = SNPs,
  genetic_variants_knockoff = snps_knockoff,
  additional_covariates = NULL,
  Z = pcs, 
  y = y,
  n_folds = 5,
  standardize = TRUE,
  FDR_rate = 0.2
)
selected_snp_names <- colnames(SNPs)[apply(results$selection_matrix, 2, any)]
colnames(results$scaled_selection_matrix) <- colnames(SNPs)
all_W_heterogeneity_plot <- plot_heatmap(
  scaled_selection_matrix = results$scaled_selection_matrix,
  selected_snps = selected_snp_names,
  ancestry = non_eur,
  ancestry_name = "non-EUR",
  category = "Selection_matrix_all",
  save_path = "~/Project/Knockoff/run_real_data_julie/marginal_test"
)


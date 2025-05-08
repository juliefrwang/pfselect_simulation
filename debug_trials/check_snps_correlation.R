library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)
library(reshape2)
library(jsonlite)

# check whether the snps are correlated 
main_dir <- "/Users/juliew/Project/Knockoff/simulation_res/seed_379_b_0.2_rep_20_true_30_hete_20_pcs_filter_MAF_diff_ans_b_clust_snps"
file_name <- glue("{main_dir}/simulated_true_data.csv")
sample_data <- read.csv(file_name)
snps <- sample_data[, grepl("chr", colnames(sample_data))]

# read beta
json_file <- glue("{main_dir}/beta_idx.json")
json_data <- fromJSON(json_file)
true_idx <- json_data$true_idx
eur_idx <- json_data$eur_idx
afr_idx <- json_data$afr_idx
homo_idx <- json_data$homo_idx

# true snps
true_snps <- snps[true_idx]

# correlation matrix 
corr_mat <- cor(true_snps)
melted_corr_mat <- melt(corr_mat)
ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2,fill=value)) + geom_tile()

ggsave(glue("snps_corr.png"), plot = last_plot(),)

# drop corr matrix that is less than 0.5
high_corr_mat <- subset(melted_corr_mat, value > 0.5) 
# heat map
ggplot(data = high_corr_mat, aes(x=Var1, y=Var2,fill=value)) + geom_tile()


beta <- numeric(ncol(snps)) 
b <- 0.2
beta[homo_idx] <- sample(c(-b, b), length(homo_idx), replace = TRUE)
beta[eur_idx] <- sample(c(-b, b), length(eur_idx), replace = TRUE)
beta[afr_idx] <- sample(c(-2*b, 2*b), length(afr_idx), replace = TRUE)

beta_matrix <- matrix(beta, nrow = nrow(snps), ncol = ncol(snps), byrow = TRUE)
#beta_matrix[, eur_idx] <- sweep(beta_matrix[, eur_idx], 1, ancestry[,'EUR'], "*")
beta_matrix[, afr_idx] <- sweep(beta_matrix[, afr_idx], 1, ancestry[,'AFR'], "*")

# y
true_y_no_noise <- rowSums(snps * beta_matrix) 

num_rows <- nrow(snps) 
fdp_list <- numeric(num_rows)
power_list <- numeric(num_rows)
power_afr_list <- numeric(num_rows)
power_homo_list <- numeric(num_rows)
num_replicate <- 20

for (i in 1:num_replicate){
  print(i)
  snps_knockoff <- generate_knockoff_data(snps)
  true_y <- true_y_no_noise + rnorm(nrow(snps))
  results <- get_importance_matrices(
    genetic_variants = snps,
    genetic_variants_knockoff = snps_knockoff,
    additional_covariates = NULL,
    Z = pcs, # or other variables like 'pcs'
    y = true_y,
    n_folds = 5,
    FDR_rate = 0.1
  )
  
  # Loop over all individuals
  for (i in 1:num_rows) {
    selected_idx_kf <- which(results$selection_matrix[i, ] == 1)
    
    # Calculate true and false positives
    num_true_pos <- sum(selected_idx_kf %in% true_idx)
    num_true_pos_afr <- sum(selected_idx_kf %in% afr_idx)
    num_true_pos_homo <- sum(selected_idx_kf %in% homo_idx)
    
    num_false_pos <- length(selected_idx_kf) - num_true_pos
    
    # Compute FDP and power
    FDP <- num_false_pos / max(1, length(selected_idx_kf))
    power <- num_true_pos / max(1, length(true_idx))
    power_afr <- num_true_pos_afr / max(1, length(afr_idx))
    power_homo <- num_true_pos_homo / max(1, length(homo_idx))
    
    # Store results
    fdp_list[i] <- fdp_list[i] + FDP
    power_list[i] <- power_list[i] + power
    power_afr_list[i] <- power_afr_list[i] + power_afr
    power_homo_list[i] <- power_homo_list[i] + power_homo
    
  }
  
}
# write.csv(cbind(fdp_list, power_list, power_eur_list, power_afr_list, power_homo_list),
#           file = glue("simulation_result_sum.csv"), row.names = FALSE)

fdp_list <- fdp_list / num_replicate
power_list <- power_list / num_replicate
power_afr_list <- power_afr_list / num_replicate
power_homo_list <- power_homo_list / num_replicate

# prepare plot
ordered_indices <- order(ancestry[,'AFR'])  # Ascending order of v
fdp_list_ordered <- fdp_list[ordered_indices]
power_list_ordered <- power_list[ordered_indices]

# create top panel: panel vs x (ordered by EUR)
power_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_list_ordered), 
             alpha = 0.7, size = 1, color = 'blue',) + 
  ylim(0,1.0) +
  theme_bw() +
  theme(panel.grid = element_blank(), ) +
  labs(x = NULL, y = "Power", color = NULL) +
  theme(legend.position = "none")

# Create bottom panel: FDR vs. x
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = fdp_list_ordered), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: EUR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices,'AFR']), 
             alpha = 0.7, size = 1, color = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Individuals ordered by ascending AFR", y = "%AFR", color = NULL)

# Combine plots
combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))

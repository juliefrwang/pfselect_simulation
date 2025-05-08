library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)

args = commandArgs(trailingOnly=TRUE)
b <- as.numeric(args[1])
seed <- as.numeric(args[2])

set.seed(seed)
num_replicate <- 100
num_true <- 30
num_affected <- 20

# check existence of working simulation folder and create if not exist
main_dir <- "/Users/juliew/Project/Knockoff"
sub_dir <- glue("simulation_res/seed_{seed}_b_{b}_rep_{num_replicate}_true_{num_true}_hete_{num_affected}_pcs_nofilter_MAF_same_b")

# check if sub directory exists 
if (file.exists(sub_dir)){
  setwd(file.path(main_dir, sub_dir))
} else {
  dir.create(file.path(main_dir, sub_dir))
  setwd(file.path(main_dir, sub_dir))
}
snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)

# Generate a 5000 x 20 matrix with random values
sample_data <- snp_data[sample(1:nrow(snp_data), 5000),] 
row.names(sample_data) <- NULL # reset the indices starting from 1.
snps <- sample_data[, grepl("chr", colnames(sample_data))]
snps <- snps[, apply(snps, 2, var) != 0]

pcs <- sample_data[,c("PC1", "PC2", "PC3", "PC4")]
# ancestry
ancestry <- as.matrix(sample_data[,c("EUR", "AFR")])

# beta
beta <- numeric(ncol(snps)) 
true_idx <- sample(1:ncol(snps), num_true)  
hete_idx <- sample(true_idx, num_affected)
homo_idx <- setdiff(true_idx, hete_idx)
eur_idx <- hete_idx[1:floor(num_affected/2)]
afr_idx <- hete_idx[(floor(num_affected/2) + 1) : num_affected]

beta[homo_idx] <- sample(c(-b, b), length(homo_idx), replace = TRUE)
beta[eur_idx] <- sample(c(-b, b), length(eur_idx), replace = TRUE)
beta[afr_idx] <- sample(c(-b, b), length(afr_idx), replace = TRUE)

beta_matrix <- matrix(beta, nrow = nrow(snps), ncol = ncol(snps), byrow = TRUE)
beta_matrix[, eur_idx] <- sweep(beta_matrix[, eur_idx], 1, ancestry[,'EUR'], "*")
beta_matrix[, afr_idx] <- sweep(beta_matrix[, afr_idx], 1, ancestry[,'AFR'], "*")
write.csv(beta_matrix, file = glue("ancestry_specified_beta.csv"), row.names = FALSE)

# y
true_y_no_noise <- rowSums(snps * beta_matrix) 

simulated_data <- cbind(snps, true_y_no_noise)
write.csv(simulated_data, file = glue("simulated_true_data.csv"), row.names = FALSE)

# Initialize result storage
num_rows <- nrow(snps)  # Number of individuals


fdp_list <- numeric(num_rows)
power_list <- numeric(num_rows)
power_eur_list <- numeric(num_rows)
power_afr_list <- numeric(num_rows)
power_homo_list <- numeric(num_rows)

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
    num_true_pos_eur <- sum(selected_idx_kf %in% eur_idx)
    num_true_pos_afr <- sum(selected_idx_kf %in% afr_idx)
    num_true_pos_homo <- sum(selected_idx_kf %in% homo_idx)
    
    num_false_pos <- length(selected_idx_kf) - num_true_pos
    
    # Compute FDP and power
    FDP <- num_false_pos / max(1, length(selected_idx_kf))
    power <- num_true_pos / max(1, length(true_idx))
    power_eur <- num_true_pos_eur / max(1, length(eur_idx))
    power_afr <- num_true_pos_afr / max(1, length(afr_idx))
    power_homo <- num_true_pos_homo / max(1, length(homo_idx))
    
    # Store results
    fdp_list[i] <- fdp_list[i] + FDP
    power_list[i] <- power_list[i] + power
    power_eur_list[i] <- power_eur_list[i] + power_eur
    power_afr_list[i] <- power_afr_list[i] + power_afr
    power_homo_list[i] <- power_homo_list[i] + power_homo
    
    if (i == 1) {
      print(FDP)
      print(power)
    }
  }
  
}

fdp_list <- fdp_list / num_replicate
power_list <- power_list / num_replicate
power_eur_list <- power_eur_list / num_replicate
power_afr_list <- power_afr_list / num_replicate
power_homo_list <- power_homo_list / num_replicate

# save output
write.csv(cbind(fdp_list, power_list, power_eur_list, power_afr_list, power_homo_list), 
          file = glue("simulation_result.csv"), row.names = FALSE)
# read output
df <- read.csv(glue("simulation_result.csv"))
fdp_list <- df$fdp_list
power_list <- df$power_list
power_eur_list <- df$power_eur_list
power_afr_list <- df$power_afr_list
power_homo_list <- df$power_homo_list

# prepare plot
ordered_indices <- order(ancestry[,'EUR'])  # Ascending order of v
fdp_list_ordered <- fdp_list[ordered_indices]
power_list_ordered <- power_list[ordered_indices]
power_eur_list_ordered <- power_eur_list[ordered_indices]
power_homo_list_ordered <- power_homo_list[ordered_indices]

ordered_indices_2 <- order(ancestry[,'AFR'])  # Ascending order of v
fdp_list_ordered_2 <- fdp_list[ordered_indices_2]
power_list_ordered_2 <- power_list[ordered_indices_2]
power_afr_list_ordered <- power_afr_list[ordered_indices_2]
power_homo_list_ordered_2 <- power_homo_list[ordered_indices_2]

# create top panel: panel vs x (ordered by EUR)
power_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_eur_list_ordered), 
             alpha = 0.7, size = 1, color = 'blue',) + 
  geom_point(aes(x = c(1:num_rows), y = power_homo_list_ordered), 
             alpha = 0.7, size = 1, color = 'orange',) + 
  ylim(0,1.0) +
  theme_bw() +
  theme(panel.grid = element_blank(), ) +
  labs(x = NULL, y = "Power", color = NULL) +
  theme(legend.position = "none")

# Create bottom panel: FDR vs. x
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = fdp_list_ordered), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey", linewidth = 0.8) + # 
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: EUR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices,'EUR']), 
             alpha = 0.7, size = 1, color = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Individuals ordered by ascending EUR", y = "%EUR", color = NULL)

# Combine plots
combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))
ggsave(glue("EUR.png"), plot = last_plot(),)


# create top panel: panel vs x (ordered by AFR) 
power_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_afr_list_ordered), 
             alpha = 0.7, size = 1, color = 'green',) + 
  geom_point(aes(x = c(1:num_rows), y = power_homo_list_ordered_2), 
             alpha = 0.7, size = 1, color = 'orange',) + 
  ylim(0,1.0) +
  theme_bw() +
  theme(panel.grid = element_blank(), ) +
  labs(x = NULL, y = "Power", color = NULL) +
  theme(legend.position = "none")

# Create middle panel: FDR vs. x
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = fdp_list_ordered_2), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey", linewidth = 0.8) + # 
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: AFR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices_2,'AFR']), 
             alpha = 0.7, size = 1, color = "grey") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Individuals ordered by ascending AFR", y = "%AFR", color = NULL)

# Combine plots
combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))
ggsave(glue("AFR.png"), plot = last_plot(),)


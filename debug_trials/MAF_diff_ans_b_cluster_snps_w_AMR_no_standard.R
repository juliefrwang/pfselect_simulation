library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)
library(jsonlite)
library(parallel)
# Set seed for reproducibility

args = commandArgs(trailingOnly=TRUE)
# b <- as.numeric(args[1])
# seed <- as.numeric(args[2])
# num_replicate <- as.numeric(args[3])
b <- 0.4
seed <- 379
num_replicate <- 100
set.seed(seed)
FDR_rate <- 0.1

num_true <- 40
num_affected <- 30

# check existence of working simulation folder and create if not exist
main_dir <- "/Users/juliew/Project/Knockoff"
sub_dir <- glue("simulation_res/seed_{seed}_b_{b}_rep_{num_replicate}_true_{num_true}_hete_{num_affected}_pcs_filter_MAF_diff_ans_b_clust_snps_w_AMR")

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
snps <- snps[, apply(snps, 2, sum) > 500]


# pcs and ancestry
pcs <- sample_data[,c("PC1", "PC2", "PC3", "PC4")]
ancestry <- as.matrix(sample_data[,c("SAS","EAS","AMR","AFR","EUR")])

# hclust snps: variables in diff clusters have  
# absolute correlation less than corr_max.
corr_mat <- cor(snps)
Sigma.distance = as.dist(1 - abs(corr_mat))
fit = hclust(Sigma.distance, method="single")
corr_max = 0.15
clusters = cutree(fit, h=1-corr_max)
cluster_rep<-match(unique(clusters),clusters) # first occurrence of each cluster in original data

# beta: select beta from cluster.index
beta <- numeric(ncol(snps)) 
true_idx <- sample(cluster_rep, num_true)  
hete_idx <- sample(true_idx, num_affected)
homo_idx <- setdiff(true_idx, hete_idx)
eur_idx <- hete_idx[1:10]
afr_idx <- hete_idx[11:20]
amr_idx <- hete_idx[21:30]


# change magnitude of the coefficient here
homo_mag <- b
eur_mag <- b
afr_mag <- b
amr_mag <- b

all_idx <- list(
  true_idx = true_idx,
  homo_idx = homo_idx,
  eur_idx = eur_idx,
  afr_idx = afr_idx,
  amr_idx = amr_idx,
  homo_mag = homo_mag,
  eur_mag = eur_mag,
  afr_mag = afr_mag,
  amr_mag = amr_mag
)
write_json(all_idx, "beta_idx.json", pretty = TRUE)


beta[homo_idx] <- sample(c(-homo_mag, homo_mag), length(homo_idx), replace = TRUE)
beta[eur_idx] <- sample(c(-eur_mag, eur_mag), length(eur_idx), replace = TRUE)
beta[afr_idx] <- sample(c(-afr_mag, afr_mag), length(afr_idx), replace = TRUE)
beta[amr_idx] <- sample(c(-amr_mag, amr_mag), length(amr_idx), replace = TRUE)

beta_matrix <- matrix(beta, nrow = nrow(snps), ncol = ncol(snps), byrow = TRUE)
beta_matrix[, eur_idx] <- sweep(beta_matrix[, eur_idx], 1, ancestry[,'EUR'], "*")
beta_matrix[, afr_idx] <- sweep(beta_matrix[, afr_idx], 1, ancestry[,'AFR'], "*")
beta_matrix[, amr_idx] <- sweep(beta_matrix[, amr_idx], 1, ancestry[,'AMR'], "*")
write.csv(beta_matrix, file = glue("ancestry_specified_beta.csv"), row.names = FALSE)

# y
true_y_no_noise <- rowSums(snps * beta_matrix) 

simulated_data <- cbind(snps,  pcs, ancestry, true_y_no_noise)
write.csv(simulated_data, file = glue("simulated_true_data.csv"), row.names = FALSE)

# Initialize result storage
num_rows <- nrow(snps)  # Number of individuals

# start R parallel computing
my_function <- function(j) {
  
  fdp_list <- numeric(num_rows)
  power_list <- numeric(num_rows)
  power_eur_list <- numeric(num_rows)
  power_afr_list <- numeric(num_rows)
  power_amr_list <- numeric(num_rows)
  power_homo_list <- numeric(num_rows)
  
  snps_knockoff <- generate_knockoff_data(snps)
  true_y <- true_y_no_noise + rnorm(nrow(snps))
  results <- get_importance_matrices(
    genetic_variants = snps,
    genetic_variants_knockoff = snps_knockoff,
    additional_covariates = NULL,
    Z = pcs, # or other variables like 'pcs'
    y = true_y,
    n_folds = 5,
    standardize = FALSE,
    FDR_rate = FDR_rate
  )
  
  json_file <- sprintf("rep_%d.json", j)
  to_save <-  list(
    coefs = as.matrix(results$coefs),
    scaled_selection_matrix = as.matrix(results$scaled_selection_matrix),
    selection_matrix = as.matrix(results$selection_matrix),
    W_statistic_matrix = as.matrix(results$W_statistic_matrix),
    q_values = as.matrix(results$q_values)
  )
  write_json(to_save, json_file, auto_unbox = TRUE)
  
  # Loop over all individuals
  for (i in 1:num_rows) {
    selected_idx_kf <- which(results$selection_matrix[i, ] == 1)
    
    # Calculate true and false positives
    num_true_pos <- sum(selected_idx_kf %in% true_idx)
    num_true_pos_eur <- sum(selected_idx_kf %in% eur_idx)
    num_true_pos_afr <- sum(selected_idx_kf %in% afr_idx)
    num_true_pos_amr <- sum(selected_idx_kf %in% amr_idx)
    num_true_pos_homo <- sum(selected_idx_kf %in% homo_idx)
    num_false_pos <- length(selected_idx_kf) - num_true_pos
    
    # Compute FDP and power
    FDP <- num_false_pos / max(1, length(selected_idx_kf))
    power <- num_true_pos / max(1, length(true_idx))
    power_eur <- num_true_pos_eur / max(1, length(eur_idx))
    power_afr <- num_true_pos_afr / max(1, length(afr_idx))
    power_amr <- num_true_pos_amr / max(1, length(amr_idx))
    power_homo <- num_true_pos_homo / max(1, length(homo_idx))
    
    # Store results
    fdp_list[i] <- fdp_list[i] + FDP
    power_list[i] <- power_list[i] + power
    power_eur_list[i] <- power_eur_list[i] + power_eur
    power_afr_list[i] <- power_afr_list[i] + power_afr
    power_amr_list[i] <- power_amr_list[i] + power_amr
    power_homo_list[i] <- power_homo_list[i] + power_homo
  }
  
  return(list(fdp_list, power_list, power_eur_list, power_afr_list, power_amr_list, power_homo_list))
}

num_cores <- detectCores() - 1
res <- mclapply(1:num_replicate, my_function, mc.cores = num_cores)

convert_to_array <- function(nested_list) {
  new_array <- simplify2array(nested_list)
  storage.mode(new_array) <- "numeric" 
  return(new_array)
}


fdp_list_all <- convert_to_array(lapply(res, `[[`, 1))
power_list_all <- convert_to_array(lapply(res, `[[`, 2))
power_eur_list_all <- convert_to_array(lapply(res, `[[`, 3))
power_afr_list_all <- convert_to_array(lapply(res, `[[`, 4))
power_amr_list_all <- convert_to_array(lapply(res, `[[`, 5))
power_homo_list_all <- convert_to_array(lapply(res, `[[`, 6))

fdp_list <- rowMeans(fdp_list_all)
power_list <- rowMeans(power_list_all)
power_eur_list <- rowMeans(power_eur_list_all)
power_afr_list <- rowMeans(power_afr_list_all)
power_amr_list <- rowMeans(power_amr_list_all)
power_homo_list <- rowMeans(power_homo_list_all)

# save output
write.csv(cbind(fdp_list, power_list, power_eur_list, power_afr_list, power_amr_list,power_homo_list), 
          file = glue("simulation_result.csv"), row.names = FALSE)

# read output
df <- read.csv(glue("simulation_result.csv"))
fdp_list <- df$fdp_list
#power_list <- df$power_list
power_eur_list <- df$power_eur_list
power_afr_list <- df$power_afr_list
power_amr_list <- df$power_amr_list
power_homo_list <- df$power_homo_list

# prepare plot
ordered_indices <- order(ancestry[,'EUR'])  # Ascending order of v
fdp_list_ordered <- fdp_list[ordered_indices]
power_eur_list_ordered <- power_eur_list[ordered_indices]
power_homo_list_ordered <- power_homo_list[ordered_indices]

ordered_indices_2 <- order(ancestry[,'AFR'])  # Ascending order of v
fdp_list_ordered_2 <- fdp_list[ordered_indices_2]
power_afr_list_ordered <- power_afr_list[ordered_indices_2]
power_homo_list_ordered_2 <- power_homo_list[ordered_indices_2]

ordered_indices_3 <- order(ancestry[,'AMR'])  # Ascending order of v
fdp_list_ordered_3 <- fdp_list[ordered_indices_3]
power_amr_list_ordered <- power_amr_list[ordered_indices_3]
power_homo_list_ordered_3 <- power_homo_list[ordered_indices_3]


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
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: EUR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices,'EUR']), 
             alpha = 0.7, size = 1, color = "grey50") +
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
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: AFR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices_2,'AFR']), 
             alpha = 0.7, size = 1, color = "grey50") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Individuals ordered by ascending AFR", y = "%AFR", color = NULL)

# Combine plots
combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))
ggsave(glue("AFR.png"), plot = last_plot(),)

# create top panel: panel vs x (ordered by AMR) 
power_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_amr_list_ordered), 
             alpha = 0.7, size = 1, color = 'cyan',) + 
  geom_point(aes(x = c(1:num_rows), y = power_homo_list_ordered_3), 
             alpha = 0.7, size = 1, color = 'orange',) + 
  ylim(0,1.0) +
  theme_bw() +
  theme(panel.grid = element_blank(), ) +
  labs(x = NULL, y = "Power", color = NULL) +
  theme(legend.position = "none")

# Create middle panel: FDR vs. x
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = fdp_list_ordered_3), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: AFR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices_3,'AMR']), 
             alpha = 0.7, size = 1, color = "grey50") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Individuals ordered by ascending AMR", y = "%AMR", color = NULL)

# Combine plots
combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))
ggsave(glue("AMR.png"), plot = last_plot(),)


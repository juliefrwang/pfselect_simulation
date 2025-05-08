# library("devtools")
# unload(package = "PFSelectcopy")
# devtools::install_local("/Users/juliew/Project/Knockoff/PFSelect_copy", force = TRUE)
setwd("/Users/juliew/Project/Knockoff/PFSelect_copy")
devtools::load_all()

library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)
library(jsonlite)
library(parallel)
library(dplyr)
# Set seed for reproducibility

args = commandArgs(trailingOnly=TRUE)

b <- 0.15
seed <- 379
num_replicate <- 100
set.seed(seed)
FDR_rate <- 0.2

# factor to consider:
sample_size <- "not_equal" #1
excl_rare <- TRUE #2
ctrl_corr <- TRUE #3
standardize <- FALSE #4
num_true <- 40 #5
num_affected <- 30
model <- "Lasso" # "Lasso" or "Grplasso" 
ctrl_ans_maf <- FALSE
# sample_size <- args[1]
# standardize <- as.logical(args[2])


# check existence of working simulation folder and create if not exist
main_dir <- "/Users/juliew/Project/Knockoff"
sub_dir <- glue("simulation_res/seed_{seed}_b_{b}_FDR_{FDR_rate}_rep_{num_replicate}_true_{num_true}_hete_{num_affected}_exclRare_{excl_rare}_ctrlCorr_{ctrl_corr}_standardize_{standardize}_sampleSize_{sample_size}_{model}_ctrlAnsMaf_{ctrl_ans_maf}_scale_pcs")

# check if sub directory exists 
if (file.exists(sub_dir)){
  setwd(file.path(main_dir, sub_dir))
} else {
  dir.create(file.path(main_dir, sub_dir))
  setwd(file.path(main_dir, sub_dir))
}
snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 

# Generate a 5000 sample
if (sample_size == "equal") {
  # Define sample sizes
  sample_sizes <- c(Admixed = 1304, EUR = 1250, AFR = 1250, AMR = 1196)
  # Create an empty list to store results
  sampled_groups <- list()
  # Process each group separately
  for (group in names(sample_sizes)) {
    group_data <- snp_data %>% filter(ans_label == group)
    size <- sample_sizes[group]
    sampled_groups[[group]] <- group_data %>% slice_sample(n = size)
  }
  # Combine all sampled groups
  sample_data <- bind_rows(sampled_groups)
} else {
  sample_data <- snp_data[sample(1:nrow(snp_data), 5000),] 
}

row.names(sample_data) <- NULL # reset the indices starting from 1.
snps <- sample_data[, grepl("chr", colnames(sample_data))]
snps <- snps[, apply(snps, 2, var) != 0]

# remove rare variants
if (excl_rare) {
  snps <- snps[, apply(snps, 2, sum) > 500]
  snps <- snps[, apply(snps, 2, sum) < 9500]
}

# pcs and ancestry
pcs <- sample_data[,c("PC1", "PC2", "PC3", "PC4")]
# scale pcs
pcs <- as.data.frame(scale(pcs))

ancestry <- as.matrix(sample_data[,c("SAS","EAS","AMR","AFR","EUR")])

# Construct snps_to_choose to select true beta: 1. ctrl_ans_maf 2.ctrl_corr
# Compute MAF within each ans_label group
if (ctrl_ans_maf) {
  for_frequency <- data.frame(snps, ans_label = sample_data$ans_label)
  for_frequency <- for_frequency %>%
    group_by(ans_label) %>%
    summarise(across(everything(), ~ pmin(mean(.x) / 2, 1 - mean(.x) / 2)))
  for_frequency <- as.data.frame(for_frequency)
  for_frequency <- for_frequency[, -1]
  for_frequency <- for_frequency[apply(for_frequency, 2, function(col) all(col > 0.1))]
  
  snps_to_choose <- snps[,colnames(for_frequency)]
  
  # 
  # var_maf <- sapply(for_frequency, var, na.rm = TRUE)
  # lowest_40_cols <- names(sort(var_maf))[1:40]
  # 
  # snps_to_choose <- snps[, lowest_40_cols]
} else {
  snps_to_choose <- snps
}


# hclust snps: variables in diff clusters have  
# absolute correlation less than corr_max.
corr_mat <- cor(snps_to_choose)
Sigma.distance = as.dist(1 - abs(corr_mat))
fit = hclust(Sigma.distance, method="single")
corr_max <- ifelse(ctrl_corr, 0.23, 1)
clusters = cutree(fit, h=1-corr_max)
cluster_rep <- match(unique(clusters),clusters) # first occurrence of each cluster in original data
cluster_rep <- colnames(snps_to_choose)[cluster_rep]

# beta: select beta from cluster.index
true_idx <- match(sample(cluster_rep, num_true) , colnames(snps))
hete_idx <- sample(true_idx, num_affected)
homo_idx <- setdiff(true_idx, hete_idx)
eur_idx <- hete_idx[1:10]
afr_idx <- hete_idx[11:20]


# change magnitude of the coefficient here
homo_mag <- eur_mag <- afr_mag <- b

# prepare for json 
all_idx <- list(
  true_idx = true_idx,
  homo_idx = homo_idx,
  eur_idx = eur_idx,
  afr_idx = afr_idx,
  homo_mag = homo_mag,
  eur_mag = eur_mag,
  afr_mag = afr_mag,
  num_snps = ncol(snps)
)
write_json(all_idx, "beta_idx.json", pretty = TRUE)

beta <- numeric(ncol(snps))
beta[homo_idx] <- sample(c(-homo_mag, homo_mag), length(homo_idx), replace = TRUE)
beta[eur_idx] <- sample(c(-eur_mag, eur_mag), length(eur_idx), replace = TRUE)
beta[afr_idx] <- sample(c(-afr_mag, afr_mag), length(afr_idx), replace = TRUE)

beta_matrix <- matrix(beta, nrow = nrow(snps), ncol = ncol(snps), byrow = TRUE)
beta_matrix[, eur_idx] <- sweep(beta_matrix[, eur_idx], 1, ancestry[,'EUR'], "*")
beta_matrix[, afr_idx] <- sweep(beta_matrix[, afr_idx], 1, ancestry[,'AFR'], "*")
write.csv(beta_matrix, file = glue("ancestry_specified_beta.csv"), row.names = FALSE)

# y
true_y_no_noise <- rowSums(snps * beta_matrix) 

simulated_data <- cbind(snps,  pcs, ancestry, true_y_no_noise)
write.csv(simulated_data, file = glue("simulated_true_data.csv"), row.names = FALSE)

# Initialize result storage
num_rows <- nrow(snps)  # Number of individuals

# create 
dir.create("Wstats")

# start R parallel computing
start_simulation <- function(j) {
  
  fdp_list <- numeric(num_rows)
  power_list <- numeric(num_rows)
  power_eur_list <- numeric(num_rows)
  power_afr_list <- numeric(num_rows)
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
    standardize = standardize,
    FDR_rate = FDR_rate,
    model = model
  )
  
  json_file <- sprintf("Wstats/rep_%d.json", j)
  to_save <-  list(
    coefs = as.vector(results$coefs),
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
    num_true_pos_homo <- sum(selected_idx_kf %in% homo_idx)
    num_false_pos <- length(selected_idx_kf) - num_true_pos
    
    # Compute FDP and power
    FDP <- num_false_pos / max(1, length(selected_idx_kf))
    power <- num_true_pos / max(1, length(eur_idx)+length(afr_idx)+length(homo_idx))
    power_eur <- num_true_pos_eur / max(1, length(eur_idx))
    power_afr <- num_true_pos_afr / max(1, length(afr_idx))
    power_homo <- num_true_pos_homo / max(1, length(homo_idx))
    
    # Store results
    fdp_list[i] <- fdp_list[i] + FDP
    power_list[i] <- power_list[i] + power
    power_eur_list[i] <- power_eur_list[i] + power_eur
    power_afr_list[i] <- power_afr_list[i] + power_afr
    power_homo_list[i] <- power_homo_list[i] + power_homo
  }
  
  return(list(fdp_list, power_list, power_eur_list, power_afr_list, power_homo_list))
}


num_cores <- detectCores() - 1

# start_simulation(1)

res <- mclapply(1:num_replicate, start_simulation, mc.cores = num_cores)

convert_to_array <- function(nested_list) {
  new_array <- simplify2array(nested_list)
  storage.mode(new_array) <- "numeric" 
  return(new_array)
}

fdp_list_all <- convert_to_array(lapply(res, `[[`, 1))
power_list_all <- convert_to_array(lapply(res, `[[`, 2))
power_eur_list_all <- convert_to_array(lapply(res, `[[`, 3))
power_afr_list_all <- convert_to_array(lapply(res, `[[`, 4))
power_homo_list_all <- convert_to_array(lapply(res, `[[`, 5))

fdp_list <- rowMeans(fdp_list_all)
power_list <- rowMeans(power_list_all)
power_eur_list <- rowMeans(power_eur_list_all)
power_afr_list <- rowMeans(power_afr_list_all)
power_homo_list <- rowMeans(power_homo_list_all)

# save output
write.csv(cbind(fdp_list, power_list, power_eur_list, power_afr_list,power_homo_list), 
          file = glue("simulation_result.csv"), row.names = FALSE)

# read output
# df <- read.csv("~/Project/Knockoff/simulation_res/seed_379_b_0.2_FDR_0.2_rep_20_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_sampleSize_not_equal_model_Grplasso_no_amr/simulation_result.csv")
df <- read.csv(glue("simulation_result.csv"))
fdp_list <- df$fdp_list
#power_list <- df$power_list
power_eur_list <- df$power_eur_list
power_afr_list <- df$power_afr_list
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
  geom_hline(yintercept = FDR_rate, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,2*FDR_rate) +
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
  geom_hline(yintercept = FDR_rate, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,2*FDR_rate) +
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


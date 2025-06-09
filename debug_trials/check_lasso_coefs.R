library(adelie)
library(glmnet)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(patchwork)
library(tidyverse)

# check LASSO whether 
sim_path <- "seed_379_b_0.23_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_sampleSize_not_equal_Lasso_ctrlAnsMaf_FALSE_only_homo_afr"
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
simulated_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
snps <- simulated_data[, grepl("chr", colnames(simulated_data))]
pcs = simulated_data[,c("PC1","PC2","PC3","PC4")]
beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx
amr_idx <- beta$amr_idx
num_snps <- ncol(snps)

# check correlation
library(corrplot)
# Extract SNP subset
homo_snps <- snps[, homo_idx]
eur_snps <- snps[, eur_idx]
afr_snps <- snps[, afr_idx]
# Set new names
colnames(homo_snps) <- paste0("homo_snp_", seq_along(homo_idx))
colnames(eur_snps) <- paste0("eur_snp_", seq_along(eur_idx))
colnames(afr_snps) <- paste0("afr_snp_", seq_along(afr_idx))
# Combine and plot
corrplot(cor(cbind(pcs, homo_snps)))
corrplot(cor(cbind(pcs, eur_snps)))
corrplot(cor(cbind(eur_snps, afr_snps)))

ggsave("corr_eur_afr.png", corrplot(cor(cbind(eur_snps, afr_snps))))
ggsave("corr_homo.png", corrplot(cor(cbind(pcs, homo_snps))))

corrplot(cor(cbind(pcs, simulated_data$true_y_no_noise)))

# get ancestry label
set.seed(379)
snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 

# Generate a 5000 sample
sample_data <- snp_data[sample(1:nrow(snp_data), 5000),] 
sample_data["non-EUR"] <- 1 - sample_data["EUR"]
# pcs and ancestry
ancestry <- as.matrix(sample_data[,c("SAS","EAS","AMR","AFR","EUR", "non-EUR")])
ordered_indices <- order(ancestry[,"non-EUR"])
ans <- sample_data[,"ans_label"]

# Check mean-W statistics over runs
mean_W <- as.matrix(read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/mean_matrix.csv")))
mean_selection <- as.matrix(read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/mean_selection_matrix.csv")))

# variance_matrix <- as.matrix(read.csv("variance_matrix.csv"))
res1.lasso <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/rep_20.json"))
W <- res1.lasso$W_statistic_matrix
  
# Step 1: Convert to long format
W_df <- as.data.frame(W)
W_df$ans <- ans
W_df$ind <- seq_len(nrow(W))

W_long <- W_df %>%
  pivot_longer(cols = -c(ans, ind), names_to = "SNP", values_to = "W")

# Step 2: Compute mean and SD per SNP per population
W_summary <- W_long %>%
  group_by(ans, SNP) %>%
  summarise(
    mean_W = mean(W),
    sd_W = sd(W),
    .groups = "drop"
  )

# Step 3: For each population, sort by abs(mean_W) and assign rank
W_summary <- W_summary %>%
  group_by(ans) %>%
  arrange(abs(mean_W), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Step 4: Lollipop plot with error bars (mean ± SD)
ggplot(W_summary, aes(x = rank, y = mean_W)) +
  geom_segment(aes(x = rank, xend = rank, y = 0, yend = mean_W), color = "lightblue") +
  geom_errorbar(aes(ymin = mean_W - sd_W, ymax = mean_W + sd_W), width = 0.2, color = "darkred") +
  geom_point(size = 1.5, color = "steelblue") +
  facet_wrap(~ans, scales = "free_x") +
  coord_cartesian(xlim = c(50, max(W_summary$rank)))+
  labs(
    title = "Mean ± SD of W-statistics by SNP (Ordered by |Mean W|)",
    x = "SNP Rank (by |Mean W|)",
    y = "Mean W-statistic"
  ) +
  theme_bw()
ggsave("W_rep_20.png", last_plot(), width = 8, height = 5)


# plot 20 replicates of mean W for AMR population
amr_subset <- which(ans == "AMR")
eur_subset <- which(ans == "EUR")
afr_subset <- which(ans == "AFR")
set.seed(2)
sample_replictes <- sample(1:100, 16)

which.subset <- afr_subset #TODO
W_long <- map_dfr(
  sample_replictes,
  function(rep_id) {
    # Load JSON file
    res <- fromJSON(glue("Wstats/rep_{rep_id}.json"))
    W <- res$W_statistic_matrix  # n x p
    sel <- res$selection_matrix
    
    # Keep only AMR samples
    W_subset <- W[which.subset,]
    sel_subset <- sel[which.subset,]
    
    # Compute mean W across AMR samples
    mean_W <- colMeans(W_subset)
    sd_W <- apply(W_subset, 2, sd)
    
    # compute thresold
    selected_flag <- apply(sel_subset, 2, function(x) any(x == 1))
    selected_prop <- colMeans(sel_subset)
    
    tibble(
      snp = paste0("SNP_", 1:ncol(W)),
      mean_W = mean_W,
      sd_W = sd_W,
      replicate = rep_id,
      snp_index = 1:ncol(W),
      highlight = snp_index %in% homo_idx,
      selected = selected_flag,
      selected_prop = selected_prop
    )
  }
)


W_long <- W_long %>%
  group_by(replicate) %>%
  arrange(abs(mean_W), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()
W_long <- W_long %>%
  mutate(
    highlight_label = ifelse(highlight, "simulated homogeneous snps", "other snps")
  )
ggplot(W_long, aes(x = rank, y = mean_W)) +
  geom_segment(aes(x = rank, xend = rank, y = 0, yend = mean_W, color = selected_prop)) +
  geom_errorbar(aes(ymin = mean_W - sd_W, ymax = mean_W + sd_W), width = 0.2, color = "darkred") +
  geom_point(aes(fill = highlight_label), shape = 21, color = "darkgrey", size = 1.6) +  # shape 21 supports fill and color
  scale_color_gradient(low = "lightblue", high = "orangered") +
  scale_fill_manual(values = c("simulated homogeneous snps" = "red", "other snps" = "white")) +
  facet_wrap(~replicate, scales = "free_x") +
  coord_cartesian(xlim = c(50, max(W_long$rank)))+
  labs(
    title = "Mean W vs SNP Rank across 16 Replicates (AFR population)",
    x = "SNP Rank (by |Mean W|)",
    y = "Mean W-statistic"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("W_afr_2.png", last_plot(), width = 12, height = 10)

# check the coefs
num_snps <- 124
for (i in seq(10:40)) {
  res <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/rep_{i}.json"))
  coefs <- unlist(res$coefs)
  x = 89
  print(coefs[5 + seq(x, x + 9 * num_snps, by = num_snps)])
}

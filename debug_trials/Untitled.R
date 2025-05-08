library(PFSelectcopy)
library(ggplot2)
library(patchwork)
library(glue)
library(jsonlite)
library(parallel)

set.seed(379)

snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)

# Generate a 5000 x 20 matrix with random values
sample_data <- snp_data[sample(1:nrow(snp_data), 5000),] 
row.names(sample_data) <- NULL # reset the indices starting from 1.
snps <- sample_data[, grepl("chr", colnames(sample_data))]
snps <- snps[, apply(snps, 2, var) != 0]
snps <- snps[, apply(snps, 2, sum) > 500]
pcs <- sample_data[,c("PC1", "PC2", "PC3", "PC4")]
ancestry <- as.matrix(sample_data[,c("EUR", "AFR")])
ordered_indices <- order(ancestry[,'EUR']) 
ordered_sample_data <- sample_data[ordered_indices]

corr_mat <- cor(snps)
melted_corr_mat <- melt(corr_mat)
ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2,fill=value)) + geom_tile()

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
eur_idx <- hete_idx[1:floor(num_affected/2)]
afr_idx <- hete_idx[(floor(num_affected/2) + 1) : num_affected]

true_snps <- snps[,true_idx]
corr_mat <- cor(true_snps)
melted_corr_mat <- melt(corr_mat)
ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2,fill=value)) + geom_tile()


# read in results with no hete effects
simulation_result <- read.csv("/Users/juliew/Project/Knockoff/simulation_res/seed_379_b_0.2_rep_25_true_filter_MAF_clust_snps_no_hete_effects/simulation_result.csv")

par(mfrow = c(2,2))
hist(simulation_result$power_list)
hist(simulation_result$power_homo_list)
hist(simulation_result$power_eur_list)
hist(simulation_result$power_afr_list)

apply(snps[,eur_idx], 2, var)
apply(snps[,homo_idx], 2, var)
apply(snps[,afr_idx], 2, var)


# plot
main_dir <- "/Users/juliew/Project/Knockoff"
sub_dir <- glue("simulation_res/seed_{seed}_b_{b}_rep_{num_replicate}_true_filter_MAF_clust_snps_all_eur_effects")
setwd(file.path(main_dir, sub_dir))
df.eur <- read.csv(glue("simulation_result.csv"))
sub_dir <- glue("simulation_res/seed_{seed}_b_{b}_rep_{num_replicate}_true_filter_MAF_clust_snps_no_hete_effects")
setwd(file.path(main_dir, sub_dir))
df.homo <- read.csv(glue("simulation_result.csv"))


plot_power_fdr <- function(fdp_list, power_list) {
  ordered_indices <- order(ancestry[,'EUR'])  # Ascending order of v
  fdp_list_ordered <- fdp_list[ordered_indices]
  power_list_ordered <- power_list[ordered_indices]
  num_rows <- length(power_list_ordered)
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
    geom_point(aes(x = c(1:num_rows), y = ancestry[ordered_indices,'EUR']), 
               alpha = 0.7, size = 1, color = "grey50") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(x = "Individuals ordered by ascending EUR", y = "%EUR", color = NULL)
  
  # Combine plots
  combined_plot <- power_plot / fdr_plot / ans_plot + plot_layout(heights = c(2, 1, 1))
  return(combined_plot) 
}

plot_eur <- plot_power_fdr(df.eur$fdp_list, df.eur$power_list) + plot_annotation(title = "All true betas interact with %EUR")
plot_homo <- plot_power_fdr(df.homo$fdp_list, df.homo$power_list) + plot_annotation(title = "All true betas stay homogenous")
combined_plots <- plot_eur | plot_homo
ggsave(glue("combined_plot.png"), plot = last_plot(),)

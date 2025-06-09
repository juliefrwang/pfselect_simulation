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

setwd("~/Project/Knockoff/run_real_data_julie") #TOCHANGE

# Load SNP data and extract snps columns
snp_filepath <- "../Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
excl_rare <- TRUE
ctrl_corr <- TRUE
which.path <- glue("run_standardize_FALSE_exclRare_TRUE_ctrlCorr_TRUE_0.9_fdr_0.2_model_Lasso_seed_400")
W_matrix <- read.csv(glue("~/Project/Knockoff/run_real_data_julie/{which.path}/W_statistic_matrix.csv"))
selection_matrix <- read.csv(glue("~/Project/Knockoff/run_real_data_julie/{which.path}/selection_matrix.csv"))


# ancestry label
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR> 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed"))))))

# define base theme
base_theme <- theme_minimal(base_size = 8) +  # Set base size for text
  theme(
    plot.title = element_text(size = 8, face = "bold"),   # Title text size
    axis.title = element_text(size = 8),                 # Axis title text size
    axis.text = element_text(size = 6),                  # Axis text size
    legend.title = element_text(size = 6),               # Legend title text size
    legend.text = element_text(size = 6)                 # Legend text size
  )

#
for_plot <- data.frame(
  PC1 <- snp_data$PC1,
  PC2 <- snp_data$PC2,
  PC3 <- snp_data$PC3,
  PC4 <- snp_data$PC4,
  ans_label <- snp_data$ans_label
)
# 
# for_plot_scaled <- data.frame(
#   PC1 <- scale(snp_data$PC1),
#   PC2 <- scale(snp_data$PC2),
#   PC3 <- scale(snp_data$PC3),
#   PC4 <- scale(snp_data$PC4),
#   ans_label <- snp_data$ans_label
# )

set2_colors <- brewer.pal(n = 8, name = "Set2")
ancestry_colors <- c("SAS" = set2_colors[1],   
                     "EAS" = set2_colors[2],  
                     "AMR" = set2_colors[3],  
                     "AFR" = set2_colors[4],  
                     "EUR" = set2_colors[5],  
                     "Admixed" = set2_colors[8]) 

# Map the colors to the ancestry labels

# pc1 and pc2
p1 <- ggplot(for_plot, aes(x = PC1, y = PC2, color = ans_label)) +
    geom_point(alpha = .7) +
    labs(title = "ADSP Cohort",
         x = "PC1",
         y = "PC2",
         color = "Ancestry label") +
    base_theme +
    scale_colour_manual(values = ancestry_colors) +
    coord_fixed(ratio = 1)

# pc3 and pc4
# p2 <- ggplot(for_plot_scaled, aes(x = PC3, y = PC4, color = ans_label)) +
#   geom_point(alpha = .7) +
#   labs(title = "ADSP Cohort PC3 vs PC4",
#        x = "PC3",
#        y = "PC4",
#        color = "Ancestry label") +
#   base_theme +
#   scale_colour_manual(values = ancestry_colors) +
#   coord_fixed(ratio = 1) 
# 
# ggsave("~/Project/Knockoff/run_real_data_julie/figs/pcs_1.png", p1)
# ggsave("~/Project/Knockoff/run_real_data_julie/figs/pcs_2.png", p2)

plot.snp.pca <- function (which.snp, snp_data, W_matrix, base_theme, excl_rare, ctrl_corr){
  # Identify which SNP column index
  x_snp <- colnames(snp_data)[grepl(which.snp, colnames(snp_data))]
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
  x_snp_idx <- which(colnames(snps) == x_snp)

  scaled_W_x <- W_matrix[, x_snp_idx] / max(W_matrix[, x_snp_idx])
  
  for_plot <- data.frame(
    PC1 <- snp_data$PC1,
    PC2 <- snp_data$PC2,
    ans_label <- snp_data$ans_label,
    scaled_W_x <- scaled_W_x
  )
  
  # for_plot_scaled <- data.frame(
  #   PC1 <- scale(snp_data$PC1),
  #   PC2 <- scale(snp_data$PC2),
  #   PC3 <- scale(snp_data$PC3),
  #   PC4 <- scale(snp_data$PC4),
  #   ans_label <- snp_data$ans_label,
  #   scaled_W_x <- scaled_W_x
  # )
  colors <- RColorBrewer::brewer.pal(11, "RdYlBu")
  x.plot <- ggplot(for_plot, aes(x = PC1, y = PC2, color = scaled_W_x)) +
    geom_point(alpha = .7) +
    labs(title = glue("{x_snp}"),
         x = "PC1",
         y = "PC2",
         color = "Scaled positive\nW-statistics") +
    base_theme +
    scale_color_gradientn(colors = colors[6:11], limits = c(0,1)) +
    coord_fixed(ratio = 1)
  return(x.plot)
}


p2 <- plot.snp.pca("APOE", snp_data, W_matrix, base_theme, excl_rare, ctrl_corr)
p3 <- plot.snp.pca("chr11.CELF1.rs10838725", snp_data, W_matrix, base_theme, excl_rare, ctrl_corr)
p4 <- plot.snp.pca("chr6.TREML2.rs60755019", snp_data, W_matrix, base_theme, excl_rare, ctrl_corr)
p5 <- plot.snp.pca("chr11.SORL1.rs11218343", snp_data, W_matrix, base_theme, excl_rare, ctrl_corr)
p6 <- plot.snp.pca("chr4.RHOH.rs2245466", snp_data, W_matrix, base_theme, excl_rare, ctrl_corr)

# combined_plot <- (p1 + p2) / (plot_spacer()+plot_spacer()) / (p3 + p4) + plot_layout(widths = c(1, -0.1, 1), heights = c(1, -0.1, 1))
# combined_plot <- p1+p2+p3+p4+plot_layout(ncol = 2)
combined_plot <- p1 + p2 + p3 + p4 + p5 + p6  #+plot_layout(widths= c(1, 1, 1))
ggsave(glue("{which.path}/fig4_pca_scaled.png"), combined_plot, width = 9, height = 4)

p1 + p2 + p3 + p4 + p5 + p6

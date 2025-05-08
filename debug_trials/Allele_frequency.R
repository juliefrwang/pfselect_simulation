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

# check existence of working simulation folder and create if not exist
main_dir <- "/Users/juliew/Project/Knockoff/"
sub_dir <- glue("simulation_res/sample_fix/")

# check if sub directory exists 
if (file.exists(sub_dir)){
  setwd(file.path(main_dir, sub_dir))
} else {
  dir.create(file.path(main_dir, sub_dir))
  setwd(file.path(main_dir, sub_dir))
}

# 1. allel frequency in different populations
# 1.1 read in simulated date and simulation results
sim_path <- "seed_379_b_0.4_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_w_AMR_sampleSize_equal"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}"))

for_plot <- data.frame(
  PC1 = simulated_true_data$PC1,
  PC2 = simulated_true_data$PC2,
  PC3 = simulated_true_data$PC3,
  PC4 = simulated_true_data$PC4,
  EUR = simulated_true_data$EUR,
  AFR = simulated_true_data$AFR,
  AMR = simulated_true_data$AMR,
  SAS = simulated_true_data$SAS,
  EAS = simulated_true_data$EAS,
  power_list = simulation_results$power_list,
  power_homo_list = simulation_results$power_homo_list,
  power_eur_list = simulation_results$power_eur_list, 
  power_afr_list = simulation_results$power_afr_list,
  power_amr_list = simulation_results$power_amr_list
)

# 1.2 allel frequency
beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx
amr_idx <- beta$amr_idx

for_frequency <- simulated_true_data[,c(homo_idx, eur_idx, afr_idx,amr_idx)]
for_plot <- for_plot %>%
  mutate(ans_label = ifelse(EUR> 0.75, "EUR",
                    ifelse(AFR > 0.75, "AFR",
                    ifelse(AMR > 0.75, "AMR",
                    ifelse(EAS > 0.75, "EAS",
                    ifelse(SAS > 0.75, "SAS", "Admixed"))))))

for_frequency$ans_label <- for_plot$ans_label

metric <- "MAF"
# Define function for MAF and Variance
compute_stat <- function(x) {
  if (metric == "MAF") {
    return(pmin(mean(x) / 2, 1 - mean(x) / 2))
  } else if (metric == "Var") {
    return(var(x))
  } else {
    stop("Invalid metric. Choose 'MAF' or 'Var'.")
  }
}
if (metric == "MAF") {
  scale_limit <- c(0, 0.5)
} else if (metric == "Var") {
  scale_limit <- c(0, 1)
}
# Compute MAF within each ans_label group
result <- for_frequency %>%
  group_by(ans_label) %>%
  summarise(across(everything(), compute_stat))
result <- as.data.frame(result)

# Compute MAF across the entire sample (no grouping)
overall_maf <- for_frequency %>%
  summarise(across(where(is.numeric), compute_stat)) %>%
  mutate(ans_label = "Overall")  # Add a label for overall MAF

# Bind the overall MAF row to the original result
result <- bind_rows(result, overall_maf)

# Convert to long format for heatmap plotting
result_long <- result %>%
  pivot_longer(-ans_label, names_to = "SNPs", values_to = "metric")
result_long$SNPs <- factor(result_long$SNPs, levels = colnames(result)[-1])  

x_axis_color <- c(rep("red", 10), rep("blue", 10), rep("purple", 10), rep("darkgreen", 10))

ggplot(result_long, aes(x = SNPs, y = ans_label, fill = metric)) +
  geom_tile() +
  geom_text(aes(label = round(metric, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient(low = "white", high = "blue", limits=scale_limit) +
  coord_fixed() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, color = x_axis_color)) +
  labs(title = glue("Heatmap of Allele {metric} by Ancestry Label"), x = "NON-NULL SNPs", y = "Ancestry")

ggsave(glue("Allele_{metric}.png"), plot = last_plot(), width = 15, height = 7)

# # 1.3 snps frequencies heat map on PC plots
# for_frequency$PC1 <- simulated_true_data$PC1
# for_frequency$PC2 <- simulated_true_data$PC2
# 
# plot_pca_by_snp <- function(for_frequency, snp_index, output_dir = "pca_by_snps") {
#   if (!dir.exists(output_dir)) dir.create(output_dir)  # Create output directory if not exist
#   col_name <- colnames(for_frequency)[snp_index]
#   ggplot(for_frequency, aes(x = PC1, y = PC2, color = .data[[col_name]])) +
#     geom_point() +
#     scale_color_gradient(low = "tomato", high = "blue", limits=c(0,2)) +  
#     theme_minimal() 
#   ggsave(filename = file.path(output_dir, paste0(i, "_", col_name, "_PCA.png")), plot = last_plot())
# }
# 
# for (i in 1:40) {
#   plot_pca_by_snp(for_frequency, i)
# }

# 2. PCA plot color based on power
power_pca <- function(for_plot, power_col) {
  plot1 <- ggplot(for_plot, aes(x = PC1, y = PC2, color = .data[[power_col]])) +
    geom_point() +
    scale_color_gradient(low = "red", high = "cadetblue1", limits=c(0,1)) +  
    theme_minimal() 
  
  plot2 <- ggplot(for_plot, aes(x = PC3, y = PC4, color = .data[[power_col]])) +
    geom_point() +
    scale_color_gradient(low = "red", high = "cadetblue1", limits=c(0,1)) +  
    theme_minimal() 
  
  combined <- (plot1 | plot2) 
  
  ggsave(glue("{power_col}-PCs.png"), 
         plot = combined, width = 12, height = 5)
}
power_pca(for_plot, "power_homo_list")
power_pca(for_plot, "power_eur_list")
power_pca(for_plot, "power_afr_list")
power_pca(for_plot, "power_amr_list")


# 3. PCA plot color based on ancestry class
for_plot <- for_plot %>%
  mutate(ans_label = ifelse(EUR> 0.75, "EUR",
                      ifelse(AFR > 0.75, "AFR",
                      ifelse(AMR > 0.75, "AMR",
                      ifelse(EAS > 0.75, "EAS",
                      ifelse(SAS > 0.75, "SAS", "Admixed"))))))

# ans.pca1 <- ggplot(for_plot, aes(x = PC1, y = PC2, color = ans_label)) +
#   geom_point() +
#   scale_colour_manual(values = c("lightgrey","blue","tomato","orchid")) +
#   theme_minimal() 
# ans.pca2 <- ggplot(for_plot, aes(x = PC3, y = PC4, color = ans_label)) +
#   geom_point() +
#   scale_colour_manual(values = c("lightgrey","blue","tomato","orchid")) +
#   theme_minimal() 
# ans.pca.combined <- (ans.pca1 | ans.pca2)
# ggsave(glue("ancestry_PCA.png"), 
#        plot = ans.pca.combined, width = 12, height = 5)
# 
# ans.single <- ggplot(for_plot, aes(x = PC1, y = PC2, color = EUR)) +
#   geom_point() +
#   scale_color_gradient(low = "white", high = "tomato", limits=c(0,1)) +  
#   theme_minimal() 
# ggsave(glue("EUR_PCA.png"), 
#        plot = ans.single, )

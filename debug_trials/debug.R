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

# 1. read in simulated date and simulation results
sim_path <- "seed_379_b_0.4_rep_100_true_40_hete_30_pcs_filter_MAF_diff_ans_b_clust_snps_w_AMR"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))

power_homo_list <- simulation_results$power_homo_list
ordered_indices <- order(simulated_true_data[,'EUR'])  # Ascending order of v

df <- data.frame(
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
  power_afr_list = simulation_results$power_afr_list
)

filtered_df <- df %>% filter(EUR == 0, power_homo_list < 0.7)

df_long <- filtered_df %>% pivot_longer(cols = c(AFR, AMR), names_to = "Population", values_to = "Value")

# Plot
ggplot(df_long, aes(x = Value, fill = Population)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue")) +
  labs(title = "Overlapping Histograms of AFR and EUR", x = "Value", y = "Count")
setwd("~/Project/Knockoff/simulation_res/seed_379_b_0.4_rep_100_true_40_hete_30_pcs_filter_MAF_diff_ans_b_clust_snps_w_AMR")
ggsave("eur0-small-power.png", plot=last_plot())

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

sim_path <- "seed_379_b_0.3_rep_100_true_50_hete_30_no_standaardize_w_AMR_sample_size_control"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}"))

# number of samples
num_rows <- nrow(simulated_true_data)

# prepare plot
ancestry_to_plot <- "HOMO" # TOCHANGE
order_by <- "EUR" # TOCHANGE
rds.file <- glue("{ancestry_to_plot}_order_by_{order_by}_beta_W_plot.rds")

# read W statistic plot
w_plot <- readRDS(glue("~/Project/Knockoff/simulation_res/{sim_path}/W_scatter/{rds.file}"))
w_plot_modified <- w_plot + xlab(NULL)

# read power
which.power <- if (ancestry_to_plot == "AMR") {
  simulation_results$power_amr_list
} else if (ancestry_to_plot == "EUR") {
  simulation_results$power_eur_list
} else if (ancestry_to_plot == "AFR") {
  simulation_results$power_afr_list
} else if (ancestry_to_plot == "HOMO") {
  simulation_results$power_homo_list
} else {
  warning("No result.")
}

# read fdp
fdp_list <- simulation_results$fdp_list

# apply order
ordered_indices <- order(simulated_true_data[, order_by])  
power_list_ordered <- which.power[ordered_indices]
fdp_list_ordered <- fdp_list[ordered_indices]

# create top panel: panel vs x (ordered by EUR)
power_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_list_ordered), 
             alpha = 0.7, size = 1, color = 'orange') + 
  ylim(0,1.0) +
  scale_x_continuous(expand = c(0.02, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Power", color = NULL) 

# Create middle panel: FDR vs. x
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = fdp_list_ordered), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  ylim(0,0.2) +
  scale_x_continuous(expand = c(0.02, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "FDR") 

# Create bottom panel: EUR vs. x 
ans_plot <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = simulated_true_data[ordered_indices,order_by]), 
             alpha = 0.7, size = 1, color = "grey50") +
  scale_x_continuous(expand = c(0.02, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = glue("Individuals ordered by ascending {order_by}"), y = (glue("%{order_by}")), color = NULL) 


# Combine plots
combined_plot <- w_plot_modified / power_plot / fdr_plot / ans_plot + plot_layout(heights = c(9,2,1,1))
ggsave(glue("{ancestry_to_plot}_beta_all.png"), plot = combined_plot, width = 12, height = 16)

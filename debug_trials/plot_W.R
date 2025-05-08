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
sim_path <- "seed_379_b_0.3_rep_100_true_50_hete_30_no_standaardize_w_AMR_sample_size_control"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}"))

file_names <- sprintf("Wstats/rep_%d.json", 1:100)

# Load all matrices
matrices <- lapply(file_names, function(file) {
  fromJSON(file)$W_statistic_matrix # Convert JSON to matrix
})

# Stack matrices into an array (5000, 126, 100)
W_array_data <- array(unlist(matrices), dim = c(dim(matrices[[1]])[1], dim(matrices[[1]])[2], length(matrices)))

# read in beta index
beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
# snps names
snps_names <- colnames(simulated_true_data)


# test 
W_test_array_data <- W_array_data[1:50, 1:12, 1:5]

# plot scatter plot of W
W.scatter.plot.df <- function(array_data, snp_idx, csv_name) {
  selected_array <- array_data[, snp_idx, ]
  df_list <- list()
  for (i in 1:length(snp_idx)) {
    df_list[[i]] <- data.frame(
          Sample = rep(1:dim(selected_array)[1], times = dim(selected_array)[3]),  # Repeat sample index for each column
          SNPs = rep(paste0("SNP-", i), times = dim(selected_array)[1] * dim(selected_array)[3]),
          Value = as.vector(selected_array[,i,])      # Flatten matrix column-wise
      )
  }
  final_df <- bind_rows(df_list)
  write.csv(final_df, csv_name)
  return(final_df)
}

dir.create(glue("~/Project/Knockoff/simulation_res/{sim_path}/W_scatter/"))
dir.create(glue("~/Project/Knockoff/simulation_res/{sim_path}/W_scatter/amr-order"))
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}/W_scatter/"))

plot_W <- function(W_array_data, ancestry_to_plot, order_by, simulated_true_data, beta, height=8) {
  # get snps index
  snp_idx <- if (ancestry_to_plot == "AMR") {
    beta$amr_idx
  } else if (ancestry_to_plot == "EUR") {
    beta$eur_idx
  } else if (ancestry_to_plot == "AFR") {
    beta$afr_idx
  } else if (ancestry_to_plot == "HOMO") {
    beta$homo_idx
  } else {
    c(beta$homo_idx, beta$eur_idx, beta$afr_idx, beta$amr_idx)
  }
  # prepare W data frame for plot
  df <- W.scatter.plot.df(W_array_data, snp_idx, glue("W_{ancestry_to_plot}_for_plot.csv"))
  level <- order(simulated_true_data[, order_by])
  df$Sample <- factor(df$Sample, levels = level)
  df$SNPs <- factor(df$SNPs, levels = sprintf("SNP-%d", 1:length(snp_idx)))
  
  # compute mean 
  df_mean <- df %>%
    group_by(SNPs, Sample) %>%
    summarize(mean_Value = mean(Value), .groups = "drop")

  
  W_plot <- ggplot(df, aes(x = Sample, y = Value)) +
    geom_point(alpha = 0.5, size = 0.5, color = "blue") +
    geom_line(data = df_mean, aes(x = Sample, y = mean_Value, group = SNPs), color = "red", linewidth = 1) +  # Mean trend line
    labs(y = "W statistic", x = glue("Sample (ascending {order_by} order)")) +
    scale_y_continuous(expand = c(0.1, 0)) +
    scale_x_discrete(expand = c(0.02, 0)) +
    facet_grid(SNPs ~ ., ) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
  saveRDS(W_plot, file=glue("{ancestry_to_plot}_order_by_{order_by}_beta_W_plot.rds"))
  ggsave(glue("{ancestry_to_plot}_order_by_{order_by}.png"), plot = W_plot, width = 12, height = height)
}

plot_W(W_array_data=W_array_data, ancestry_to_plot="ALL", order_by="EUR", simulated_true_data=simulated_true_data, beta=beta, height = 20)

plot_W(W_array_data=W_array_data, ancestry_to_plot="AMR", order_by="AMR", simulated_true_data=simulated_true_data, beta=beta)
plot_W(W_array_data=W_array_data, ancestry_to_plot="EUR", order_by="EUR", simulated_true_data=simulated_true_data, beta=beta)
plot_W(W_array_data=W_array_data, ancestry_to_plot="AFR", order_by="AFR", simulated_true_data=simulated_true_data, beta=beta)
plot_W(W_array_data=W_array_data, ancestry_to_plot="HOMO", order_by="EUR", simulated_true_data=simulated_true_data, beta=beta, height=16)

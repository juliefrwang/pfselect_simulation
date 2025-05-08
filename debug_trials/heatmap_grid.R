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

sim_path <- "seed_379_b_0.4_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_w_AMR_sampleSize_equal"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))

beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx
amr_idx <- beta$amr_idx

setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats"))
file_names <- sprintf("rep_%d.json", 1:100)

# Load all matrices
matrices <- lapply(file_names, function(file) {
  fromJSON(file)$W_statistic_matrix # Convert JSON to matrix
})

# Stack matrices into an array (5000, 126, 100)
array_data <- array(unlist(matrices), dim = c(dim(matrices[[1]])[1], dim(matrices[[1]])[2], length(matrices)))

# Compute variance across repetitions
variance_matrix <- apply(array_data, c(1, 2), var)
mean_matrix <- apply(array_data, c(1, 2), mean)

# Save variance matrix to CSV
write.csv(variance_matrix, "variance_matrix.csv", row.names = FALSE)
write.csv(mean_matrix, "mean_matrix.csv", row.names = FALSE)

# Load all matrices
matrices <- lapply(file_names, function(file) {
  fromJSON(file)$selection_matrix # Convert JSON to matrix
})

# Stack matrices into an array (5000, 126, 100)
array_data <- array(unlist(matrices), dim = c(dim(matrices[[1]])[1], dim(matrices[[1]])[2], length(matrices)))

# Compute variance across repetitions
var_selection_matrix <- apply(array_data, c(1, 2), var)
mean_selection_matrix <- apply(array_data, c(1, 2), mean)

# Save variance matrix to CSV
write.csv(var_selection_matrix, "var_selection_matrix.csv", row.names = FALSE)
write.csv(mean_selection_matrix, "mean_selection_matrix.csv", row.names = FALSE)

simulated_true_data$non_EUR <- 1 - simulated_true_data$EUR
ordered_indices <- order(simulated_true_data[,"non_EUR"])

# define plot function
heatmap.plot <- function(matrix_name, matrix, index, grid_label_color, ordered_indices = ordered_indices) {
  matrix_for_plot <- matrix[,index]
  ordered_map <- matrix_for_plot[ordered_indices,]
  df_melt <- reshape2::melt(ordered_map)
  colnames(df_melt) <- c("Samples", "SNPs", matrix_name)
  df_melt$grid <- matrix_name
  
  # Plot var ord heat map
  plot_ <- ggplot(df_melt, aes(x = as.factor(Samples), y = as.factor(SNPs), fill = df_melt[,matrix_name])) +
    geom_tile() +
    scale_fill_gradient(low = "aliceblue", high = "orangered", name = matrix_name) +
    scale_x_discrete(
      breaks = c("1", "1000", "2000", "3000", "4000", "5000"),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_grid(rows = vars(grid), scales = "free", space = "free") +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = grid_label_color, color = NA),
      plot.margin = margin(0,0,0,0)
    ) +
    labs(x = NULL, y = "SNPs") 
  
  return(plot_)
}


which.idx <- homo_idx # CONFIG 1:

# 1. Three Heatmaps 
var.w.plot <- heatmap.plot("Variance W", variance_matrix, c(which.idx), "lightpink")
mean.w.plot <- heatmap.plot("Mean W", mean_matrix, c(which.idx), "lavender")
var.selection.plot <- heatmap.plot("Variance of selection", var_selection_matrix, c(which.idx), "lightpink")
mean.selection.plot <- heatmap.plot("Mean of selection", mean_selection_matrix, c(which.idx), "lavender") + 
  scale_fill_gradient(low = "aliceblue", high = "orangered", name = "Mean of selection", limits = c(0,1))

ans_plot <- ggplot() +
  geom_point(aes(x = c(1:5000), 
                 y = simulated_true_data[ordered_indices,"non_EUR"],
                 color = "non_EUR"), 
             alpha = 0.7, size = 0.6) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(1, 1000, 2000, 3000, 4000, 5000)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(values = c("non_EUR" = "darkgrey"), 
                     labels = c("non_EUR" = "non European")) +
  labs(x = glue("Samples (ordered by asceding non-EUR)"), y = "Percentage", color = "Group")

# 2. FDR Plot
# Create middle panel: FDR vs. x
fdp_list <- simulation_results$fdp_list
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:5000), y = fdp_list[ordered_indices]), 
             alpha = 0.7, size = 1, color = "purple") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(1, 1000, 2000, 3000, 4000, 5000)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.2)) +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.position = "none",) +
  labs(x = NULL, y = "FDR")

# 3. Ancestry Plot
which.ans <- "EUR" # CONFIG 2: "EUR", "AFR", "AMR" or NULL
if (!is.null(which.ans)) {
  ans_plot <- ans_plot +
    geom_point(aes(x = c(1:5000), 
                   y = simulated_true_data[ordered_indices, which.ans],  
                   color = which.ans),  # Color is mapped to the value of `which.ans`
               alpha = 0.7, size = 0.6)
  ans_plot <- ans_plot +
    scale_color_manual(values = c("EUR" = "blue", "non_EUR" = "darkgrey"), 
                       labels = c("EUR" = "European", "non_EUR" = "non European")) +
    labs(color = "Group") 
}



# Combine plots
combined_plot <- mean.w.plot / var.w.plot / mean.selection.plot / fdr_plot/ ans_plot + plot_annotation(tag_levels = 'a') 
ggsave(glue("../combined_heatmaps_homo_2.png"), combined_plot, width = 9, height = 8)


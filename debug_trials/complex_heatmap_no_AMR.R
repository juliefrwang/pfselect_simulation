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
library(cowplot)
library(ggtext)

standardize <-FALSE
sim_path <- "seed_379_b_0.23_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_sampleSize_not_equal_Lasso_ctrlAnsMaf_FALSE"
simulated_true_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))

beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx


setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats"))

# generate W and selection heatmap matrix 
# read in files
file_names <- sprintf("rep_%d.json", 1:20)
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


########## Define Heatmap Function ###########
### 1. W statistics or selection matrix heatmap ###
heatmap.plot.1 <- function(matrix_name, matrix, which.idx, ordered_indices) {
  matrix_for_plot <- matrix[,which.idx]
  
  # order tha idx by mu/sd
  snp_std_mu <- apply(matrix_for_plot, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  matrix_for_plot <- matrix_for_plot[, order(snp_std_mu)]
  # end
  
  ordered_map <- matrix_for_plot[ordered_indices,]
  df_melt <- reshape2::melt(ordered_map)
  colnames(df_melt) <- c("Samples", "SNPs", matrix_name)
  
  # custom_limits <- if (grepl("selection", matrix_name)) c(0, 1) else NULL
  if (grepl("selection", matrix_name)) {
    custom_limits <- c(0, 1)
  } else {
    custom_limits <- c(0, max(df_melt[[matrix_name]], na.rm = TRUE))
  }
  legend_name <- if (grepl("selection", matrix_name)) "Mean of Selection" else matrix_name
  
  # Plot var ord heat map
  plot_ <- ggplot(df_melt, aes(x = as.factor(Samples), y = as.factor(SNPs), fill = df_melt[,matrix_name])) +
    geom_tile() +
    scale_fill_gradient(low = "aliceblue", high = "orangered", name = legend_name, limits = custom_limits) +
    scale_x_discrete(
      breaks = c("1", "1000", "2000", "3000", "4000", "5000"),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic(base_size = 8) +
    guides(fill = guide_colorbar(keywidth = 0.3, keyheight = 2)) + # to the right of the plot
    theme(
      legend.title = element_text(size = 7), 
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 8),
      plot.margin = margin(6,6,0,0),
    ) +
    labs(x = NULL, y = "SNPs", title = matrix_name) 
  
  return(plot_)
}

heatmap.plot <- function(matrix_name, matrix, which.idx, ordered_indices) {
  matrix_for_plot <- matrix[,which.idx]
  
  # order tha idx by mu/sd
  snp_std_mu <- apply(matrix_for_plot, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  matrix_for_plot <- matrix_for_plot[, order(snp_std_mu)]
  # end
  
  ordered_map <- matrix_for_plot[ordered_indices,]
  df_melt <- reshape2::melt(ordered_map)
  colnames(df_melt) <- c("Samples", "SNPs", matrix_name)
  
  # custom_limits <- if (grepl("selection", matrix_name)) c(0, 1) else NULL
  if (grepl("selection", matrix_name)) {
    custom_limits <- c(0, 1)
  } else {
    custom_limits <- c(0, max(df_melt[[matrix_name]], na.rm = TRUE))
  }
  legend_name <- if (grepl("selection", matrix_name)) "Mean of Selection" else matrix_name
  
  # Plot var ord heat map
  plot_ <- ggplot(df_melt, aes(x = as.factor(Samples), y = as.factor(SNPs), fill = df_melt[,matrix_name])) +
    geom_tile() +
    scale_fill_gradient(low = "aliceblue", high = "orangered", name = legend_name, limits = custom_limits) +
    scale_x_discrete(
      breaks = c("1", "1000", "2000", "3000", "4000", "5000"),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic(base_size = 8) +
    guides(fill = guide_colorbar(keywidth = 5, keyheight = 0.3)) + # to the bottom of the plot
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.position = "bottom",
      plot.title = element_text(size = 8, face = "bold"),
      plot.margin = margin(4,4,0,0)
    ) +
    labs(x = NULL, y = "SNPs", title = matrix_name) 
  
  return(plot_)
}
### 2.Allele frequency heatmap ###
allele.plot <- function(simulated_true_data, which.idx) {
  for_frequency <- simulated_true_data[,c(which.idx)]
  for_plot <- simulated_true_data %>%
    mutate(ans_label = ifelse(EUR> 0.75, "EUR",
                       ifelse(AFR > 0.75, "AFR",
                       ifelse(AMR > 0.75, "AMR",
                       ifelse(EAS > 0.75, "EAS",
                       ifelse(SAS > 0.75, "SAS", "Admixed"))))))
  
  for_frequency$ans_label <- for_plot$ans_label
  
  # Compute MAF within each ans_label group
  result <- for_frequency %>%
    group_by(ans_label) %>%
    summarise(across(everything(), ~ pmin(mean(.x) / 2, 1 - mean(.x) / 2)))
  result <- as.data.frame(result)
  
  # Compute MAF across the entire sample (no grouping)
  overall_maf <- for_frequency %>%
    summarise(across(where(is.numeric), ~ pmin(mean(.x) / 2, 1 - mean(.x) / 2))) %>%
    mutate(ans_label = "Overall")  # Add a label for overall MAF
  
  # Bind the overall MAF row to the original result
  result <- bind_rows(result, overall_maf)
  
  # Convert to long format for heatmap plotting
  result_long <- result %>%
    pivot_longer(-ans_label, names_to = "SNPs", values_to = "MAF")
  result_long$SNPs <- factor(result_long$SNPs, levels = colnames(result)[-1])  
  
  allele_plot <- ggplot(result_long, aes(x = ans_label, y = SNPs, fill = MAF)) +
    geom_tile() +
    geom_text(aes(label = round(MAF, 2)), color = "black", size = 2) +  # Add text labels
    scale_fill_gradient(low = "white", high = "blue", limits=c(0, 0.5)) +
    theme_minimal(base_size = 8) +
    guides(fill = guide_colorbar(keywidth = 0.3, keyheight = 2)) + 
    theme(
      # axis.ticks.y = element_blank(),
      # axis.text.y = element_blank(),
      legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6),
      legend.position = "none",
      plot.title = element_text(size = 8)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "MAF")
  
  return(allele_plot)
}


mean_matrix <- as.matrix(read.csv("mean_matrix.csv"))
# variance_matrix <- as.matrix(read.csv("variance_matrix.csv"))
mean_selection_matrix <- as.matrix(read.csv("mean_selection_matrix.csv"))
# var_selection_matrix <- as.matrix(read.csv("var_selection_matrix.csv"))
colnames(mean_matrix) <- NULL
# colnames(variance_matrix) <- NULL
colnames(mean_selection_matrix) <- NULL
# colnames(var_selection_matrix) <- NULL

# define index ordered by non-european
simulated_true_data$non_EUR <- 1 - simulated_true_data$EUR
ordered_indices <- order(simulated_true_data[,"non_EUR"])

########## START  PLOTING ##############

########## 1. Appendix Heatmap + MAF ##############
heatmap_w_maf <- function(which.idx, file.ans.tag) {
  
  ### Heatmaps ###
  snp_label <- as.character(1:10) 
  if (file.ans.tag == "HOMO") {
    snp_label <- as.character(1:10)
  } else if (file.ans.tag == "EUR") {
    snp_label <- as.character(11:20)
  } else if (file.ans.tag == "AFR") {
    snp_label <- as.character(21:30)
  }
  mean.w.plot <- heatmap.plot.1("Mean W", mean_matrix, c(which.idx), ordered_indices) + scale_y_discrete(labels = snp_label)
  mean.selection.plot <- heatmap.plot.1("Mean selection", mean_selection_matrix, c(which.idx), ordered_indices) + scale_y_discrete(labels = snp_label)
      
  ### Allele Frequency Plot ###
  allele.freq.plot <- allele.plot(simulated_true_data, which.idx) + scale_y_discrete(labels = snp_label)
  
  ### Combine Plot ###
  mean.w.allele.plot <- mean.w.plot + allele.freq.plot + plot_layout(widths = c(3,1), guides = "collect") 
  mean.selection.allele.plot <- mean.selection.plot + allele.freq.plot + plot_layout(widths = c(3,1), guides = "collect") 
  
  combined_plot <- mean.w.allele.plot / mean.selection.allele.plot +
    plot_annotation(tag_levels = 'a', title = glue("standardize = {standardize}, excluding rare variants, FDR = 0.2\n{file.ans.tag}") )
  
  ggsave(glue("../heatmaps_maf_{file.ans.tag}_2.png"), combined_plot, width = 7.5, height = 3.5, dpi = 300, units = "in")
  
}
heatmap_w_maf(eur_idx, "EUR")
heatmap_w_maf(afr_idx, "AFR")
heatmap_w_maf(homo_idx, "HOMO")

############## 2. All heatmaps ##############

### FDR Plot  ###
# Create middle panel: FDR vs. x
fdp_list <- simulation_results$fdp_list
FDR_thresold <- 0.2
# if (abs(0.2-mean(fdp_list)) < abs(0.1-mean(fdp_list))) {FDR_thresold = 0.2} else {FDR_thresold = 0.1} 
fdr_plot <- ggplot() +
  geom_point(aes(x = c(1:5000), y = fdp_list[ordered_indices]), 
             alpha = 0.7, size = 0.6, color = "purple") +
  geom_hline(yintercept = FDR_thresold, linetype = "dashed", color = "grey50", 
             linewidth = 0.8) + 
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(1, 1000, 2000, 3000, 4000, 5000)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,2*FDR_thresold)) +
  theme_classic(base_size = 8) +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(4,4,0,0)) +
  labs(x ="Samples (ordered by ascending non-EUR)", y = "FDR")

### Ancestry Plot  ###
set_colors <- brewer.pal(n = 8, name = "Set2")
ancestry_colors <- c("SAS" = set_colors[1],   
                     "EAS" = set_colors[2],  
                     "AMR" = set_colors[3],  
                     "AFR" = set_colors[4],  
                     "EUR" = set_colors[5],  
                     "Admixed" = set_colors[8]) 

ancestry_plot <- function(which.ans) {
  if (!is.null(which.ans) && which.ans == "EUR") {
    legend.p <- c(0.25, 0.5)
  } else {
    legend.p <- c(0.25, 0.8)
    }
  if (!is.null(which.ans)) {
    ans_plot <- ggplot() +
      geom_point(aes(x = c(1:5000), 
                     y = simulated_true_data[ordered_indices,"non_EUR"],
                     color = "non_EUR"), 
                 alpha = 0.7, size = 0.3) +
      geom_point(aes(x = c(1:5000), 
                     y = simulated_true_data[ordered_indices, which.ans],  
                     color = which.ans),  # Color is mapped to the value of `which.ans`
                 alpha = 0.5, size = 0.3) +
      scale_x_continuous(expand = c(0, 0),
                         breaks = c(1, 1000, 2000, 3000, 4000, 5000)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_classic(base_size = 8) +
      theme(
        legend.position = legend.p,
        legend.title = element_blank(), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(1, 'line'),
        legend.background = element_rect(fill=NULL, colour = NULL),
        plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(4,4,0,0)
      ) + 
      scale_color_manual(values = c("AMR" = set_colors[3], "EUR" = set_colors[5], "AFR" = set_colors[4], "non_EUR" = "red"), 
                         labels = c(which.ans = which.ans, "non_EUR" = "non-EUR")) +
      labs(x = glue("Samples \nordered by ascending non-EUR"), y = "Percent", color = "Group")
    
  } else {
    ans_plot <- ggplot() +
      geom_point(aes(x = c(1:5000), 
                     y = simulated_true_data[ordered_indices,"non_EUR"],
                     color = "non_EUR"), 
                 alpha = 0.7, size = 0.3) +
      scale_x_continuous(expand = c(0, 0),
                         breaks = c(1, 1000, 2000, 3000, 4000, 5000)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_classic(base_size = 8) +
      theme(
        legend.position = legend.p,
        legend.title = element_blank(), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(1, 'line'),
        legend.background = element_rect(fill=NULL, colour = NULL),
        plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(4,4,0,0)
      ) + 
      scale_color_manual(values = c("AMR" = set_colors[3], "EUR" = set_colors[5], "AFR" = set_colors[4], "non_EUR" = "red"), 
                         labels = c(which.ans = which.ans, "non_EUR" = "non-EUR")) +
      labs(x = glue("Samples \nordered by ascending non-EUR"), y = "Percent", color = "Group", title = "Ancestry Composition")
  }
  return(ans_plot)
}

### prepare all sub plots ###
mean.w.plot.1 <- heatmap.plot("Mean W \nHomogeneous effect size", mean_matrix, homo_idx, ordered_indices) + scale_y_discrete(labels = as.character(1:10))
mean.selection.plot.1 <- heatmap.plot("Mean selection \nHomogeneous effect size", mean_selection_matrix, homo_idx, ordered_indices) + scale_y_discrete(labels = as.character(1:10))

mean.w.plot.2 <- heatmap.plot("Mean W \nEffect size interacting with EUR%", mean_matrix, eur_idx, ordered_indices) + scale_y_discrete(labels = as.character(11:20))
mean.selection.plot.2 <- heatmap.plot("Mean selection \nEffect size interacting with EUR%", mean_selection_matrix, eur_idx, ordered_indices) + scale_y_discrete(labels = as.character(11:20))

mean.w.plot.3 <- heatmap.plot("Mean W \nEffect size interacting with AFR%", mean_matrix, afr_idx, ordered_indices) + scale_y_discrete(labels = as.character(21:30))
mean.selection.plot.3 <- heatmap.plot("Mean selection \nEffect size interacting with AFR%", mean_selection_matrix, afr_idx, ordered_indices) + scale_y_discrete(labels = as.character(21:30))


ancestry_plot.1 <- ancestry_plot(NULL) 
ancestry_plot.2 <- ancestry_plot("EUR") 
ancestry_plot.3 <- ancestry_plot("AFR") 

layout <- c(
  area(1, 1), area(1, 2), area(1, 3),
  area(2, 1), area(2, 2), area(2, 3), 
  area(3, 1), area(3, 2), area(3, 3), 
  area(4, 1, 4,3) # guide_area spans two columns
)

combined_plot_2 <- mean.w.plot.1 + mean.w.plot.2 + mean.w.plot.3 +
  mean.selection.plot.1 + mean.selection.plot.2 + mean.selection.plot.3 +
  ancestry_plot.1 + ancestry_plot.2 + ancestry_plot.3  +
  fdr_plot + plot_layout(design = layout) + plot_annotation(title = glue("standardize = {standardize}, excluding rare variants,  FDR = 0.2"), tag_levels = 'a') &  theme(plot.tag = element_text(size = 8,hjust = 5, vjust = 0))
 
ggsave("../combined_plot_test8.png", combined_plot_2, width = 8, height = 7, units = "in", dpi = 400)


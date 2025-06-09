# Load the package
library(PFSelectcopy)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Load matrices
load_matrices <- function(save_path) {
    scaled_selection_matrix <- read.csv(file.path(save_path, "scaled_selection_matrix.csv"))
    selection_matrix <- read.csv(file.path(save_path, "selection_matrix.csv"))
    W_statistic_matrix <- read.csv(file.path(save_path, "W_statistic_matrix.csv"))
    q_values <- read.csv(file.path(save_path, "q_values.csv"))
    results <- list(scaled_selection_matrix = scaled_selection_matrix,
       selection_matrix = selection_matrix,
       W_statistic_matrix = W_statistic_matrix,
       q_values = q_values)
    return(results)
}

# Categorize SNPs by standard deviation of importance scores
categorize_snps_by_std_dev <- function(scaled_selection_matrix) {
  # Compute standard deviation for each SNP
  snp_std_dev <- apply(scaled_selection_matrix, 2, sd, na.rm = TRUE)
  
  # Rank SNPs by standard deviation (ascending)
  ranked_snps <- names(sort(snp_std_dev))
  n <- length(ranked_snps)
  
  # Categorize SNPs
  snp_categories <- list(
    low_heterogeneity = ranked_snps[1:(0.35 * n)],  # Lowest 35% by standard deviation
    medium_heterogeneity = ranked_snps[(0.35 * n + 1):(0.7 * n)],  # Middle 35%
    high_heterogeneity = ranked_snps[(0.7 * n + 1):n]  # Top 30%
  )
  return(snp_categories)
}

# Print standard deviation for SNPs in each category
print_snp_std_dev <- function(snp_categories, scaled_selection_matrix) {
  for (category in names(snp_categories)) {
    cat("\n", category, "category SNPs and their standard deviations:\n")
    
    # Extract SNPs for the category
    snps <- snp_categories[[category]]
    
    # Get standard deviation for these SNPs
    snp_std_devs <- apply(scaled_selection_matrix[, snps, drop = FALSE], 2, sd, na.rm = TRUE)
    
    # Print SNPs and their standard deviations
    sorted_std_devs <- sort(snp_std_devs)
    print(sorted_std_devs)
  }
}


select_top_snps <- function(snp_category, category_name, num_snps = 10) {
  selected_snps <- tail(snp_category, num_snps) 

  cat(paste("Selected SNPs for", category_name, "heterogeneity:\n"))
  print(selected_snps)
  
  if (length(selected_snps) < num_snps) {
    warning(paste("Only", length(selected_snps), "SNPs available in this category."))
  }
  return(selected_snps)
}

plot_heatmap <- function(scaled_selection_matrix, selected_snps, ancestry, ancestry_name, category, save_path) {
  # Extract chromosome numbers for sorting
  chromosome_order <- as.numeric(gsub("chr([0-9]+)\\..*", "\\1", selected_snps))
  snp_data <- data.frame(SNP = selected_snps, Chromosome = chromosome_order)
  
  # Sort selected SNPs by chromosome for visualization
  snp_data_sorted <- snp_data[order(snp_data$Chromosome), ]
  sorted_snps <- snp_data_sorted$SNP
  
  # Filter matrix for selected SNPs
  scaled_selection_matrix_subset <- scaled_selection_matrix[, sorted_snps, drop = FALSE]
  
  # Sort individuals by "ancestry" percentage
  sorted_indices <- order(ancestry)
  sorted_scaled_matrix <- scaled_selection_matrix_subset[sorted_indices, ]
  
  # Convert matrix to long format for ggplot
  long_matrix <- melt(as.matrix(sorted_scaled_matrix))
  colnames(long_matrix) <- c("Individual", "SNP", "Importance")
  
  # Ensure Individuals are sorted by EUR for ggplot
  long_matrix$Individual <- factor(long_matrix$Individual, levels = sorted_indices)  # Align with EUR sort order
  
  # Shorten SNP names for y-axis labels (retain chr and SNP name only)
  short_snps <- gsub("(chr[0-9]+\\.[^\\.]+)\\..*", "\\1", sorted_snps)
  # Ensure SNPs are sorted by chromosome in the heatmap
  # long_matrix$SNP <- factor(long_matrix$SNP, levels = sorted_snps, labels = short_snps)
  long_matrix$SNP <- factor(long_matrix$SNP, levels = sorted_snps)
  
  # Handle cases with a single unique non-zero value
  non_zero_importance <- long_matrix$Importance[long_matrix$Importance > 0]
  if (length(unique(non_zero_importance)) == 1) {
    gradient_values <- c(0, 1)  # Default gradient for single unique value
  } else {
    gradient_values <- scales::rescale(c(0, min(non_zero_importance), max(non_zero_importance)))
  }
  
  # Define color palette and scaling logic
  colors <- RColorBrewer::brewer.pal(11, "RdYlBu")
  heatmap_plot <- ggplot(long_matrix, aes(x = Individual, y = SNP, fill = Importance)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("lightgrey", colors),
      values = gradient_values,
      na.value = "lightgrey", 
      limits = c(0, max(non_zero_importance, na.rm = TRUE))
    ) +
    labs(
      title = paste(category, "Heterogeneity SNPs"),
      x = paste("Individuals (Ranked by", ancestry_name, "%)"),
      y = "Genetic Variants (SNP)",
      fill = "Relative\nImportance"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold", margin = margin(r = 10)),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
      legend.title = element_text(size = 7), 
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )
  
  # Save plot
  file_name <- file.path(save_path, paste0(ancestry_name, "_", category, "_heatmap_sorted.png"))
  ggsave(file_name, plot = heatmap_plot, width = 7.5, height = 4, units = "in")
  print(paste(category, "heatmap saved to", file_name))
  return(heatmap_plot)
}



# main section: 
setwd("~/Project/Knockoff/run_real_data_julie")
snp_filepath <- "~/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv" # rows: samples, cols: PCs, EUR, ..., SNPs, Outcome(AD)
snp_data <- read.csv(snp_filepath)

pcs <- snp_data[, c("PC1", "PC2", "PC3", "PC4")]
eur <- snp_data$EUR
afr <- snp_data$AFR
amr <- snp_data$AMR
y <- snp_data$AD
genetic_variants <- snp_data[, grepl("chr", colnames(snp_data))] 

load_directory = "~/Project/Knockoff/run_real_data_julie/runs_standardize_TRUE_ctrlMAF_FALSE_fdr_0.1_pcs_interacting"
matrices <- load_matrices(load_directory)
colnames(matrices$scaled_selection_matrix) <- colnames(genetic_variants)
colnames(matrices$selection_matrix) <- colnames(genetic_variants)
colnames(matrices$W_statistic_matrix) <- colnames(genetic_variants)

selected_snp_names <- colnames(genetic_variants)[apply(matrices$selection_matrix, 2, any)]
print("Selected features based on the knockoff filter:")
print(selected_snp_names)

# Filter scaled_selection_matrix to include only selected SNPs
scaled_selection_matrix_filtered <- matrices$scaled_selection_matrix[, selected_snp_names, drop = FALSE]
# Categorize SNPs by non-zero counts using the filtered matrix
snp_categories <- categorize_snps_by_std_dev(scaled_selection_matrix_filtered)
print_snp_std_dev(snp_categories, scaled_selection_matrix_filtered)

low_snps <- select_top_snps(snp_categories$low_heterogeneity, "low", num_snps = 10)
medium_snps <- select_top_snps(snp_categories$medium_heterogeneity, "medium", num_snps = 10)
high_snps <- select_top_snps(snp_categories$high_heterogeneity, "high", num_snps = 10)

save_directory <- "figs" 
non_eur <- 1 - eur
# low_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = low_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Low",
#   save_path = save_directory
# )
# 
# medium_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = medium_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Medium",
#   save_path = save_directory
# )

high_heterogeneity_plot <- plot_heatmap(
  scaled_selection_matrix = matrices$scaled_selection_matrix,
  selected_snps = high_snps,
  ancestry = non_eur,
  ancestry_name = "non-EUR",
  category = "High",
  save_path = save_directory
)


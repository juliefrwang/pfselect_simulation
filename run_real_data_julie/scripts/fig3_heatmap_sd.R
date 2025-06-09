# Load the package
library(PFSelectcopy)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(glue)


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
    ranked_snps = ranked_snps,
    low_heterogeneity = ranked_snps[1:ceiling(0.33*n)],  # Lowest 35% by standard deviation
    medium_heterogeneity = ranked_snps[(ceiling(0.33*n) + 1):ceiling(0.67 * n)],  # Middle 35%
    high_heterogeneity = ranked_snps[(ceiling(0.67 * n) + 1):n]  # Top 30%
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
  # selected_snps <- tail(snp_category, num_snps) 
  selected_snps <- head(snp_category, num_snps)

  cat(paste("Selected SNPs for", category_name, "heterogeneity:\n"))
  print(selected_snps)
  
  if (length(selected_snps) < num_snps) {
    warning(paste("Only", length(selected_snps), "SNPs available in this category."))
  }
  return(selected_snps)
}

plot_heatmap <- function(scaled_selection_matrix, selected_snps, ancestry, ancestry_name, category, save_path, max_non_zero=max_non_zero, min_non_zero=min_non_zero, w, h) {
  if (grepl("ranked", category)) {
    scaled_selection_matrix_subset <- scaled_selection_matrix[, selected_snps, drop = FALSE]
    sorted_snps <- selected_snps
  } else {
    # Extract chromosome numbers for sorting
    chromosome_order <- as.numeric(gsub("chr([0-9]+)\\..*", "\\1", selected_snps))
    snp_data <- data.frame(SNP = selected_snps, Chromosome = chromosome_order)
    
    # Sort selected SNPs by chromosome for visualization
    snp_data_sorted <- snp_data[order(snp_data$Chromosome), ]
    sorted_snps <- snp_data_sorted$SNP
    
    # Filter matrix for selected SNPs
    scaled_selection_matrix_subset <- scaled_selection_matrix[, sorted_snps, drop = FALSE]
  }
  
  # Sort individuals by "ancestry" percentage
  sorted_indices <- order(ancestry)
  sorted_scaled_matrix <- scaled_selection_matrix_subset[sorted_indices, ]
  rownames(sorted_scaled_matrix) <- NULL
  
  # Convert matrix to long format for ggplot
  long_matrix <- reshape2::melt(as.matrix(sorted_scaled_matrix))
  colnames(long_matrix) <- c("Individual", "SNP", "Importance")
  
  # Ensure SNPs are sorted by chromosome in the heatmap
  long_matrix$SNP <- factor(long_matrix$SNP, levels = sorted_snps)
  
  # Handle cases with a single unique non-zero value
  non_zero_importance <- long_matrix$Importance[long_matrix$Importance > 0]
  
  if (length(unique(non_zero_importance)) == 1) {
    gradient_values <- c(0, 1)  # Default gradient for single unique value
    print(gradient_values)
  } else {
    gradient_values <- scales::rescale(c(0, min_non_zero, max_non_zero))
    gradient_values <- c(0, seq(min_non_zero, max_non_zero, (max_non_zero-min_non_zero)/5))
    # print(gradient_values)
    # gradient_values <- sort(gradient_values)
    # print(gradient_values)
    # gradient_values <- c(0, seq(gradient_values[2], gradient_values[3], (gradient_values[3]-gradient_values[2])/5))
    # print(gradient_values)
  }
  
  # Define color palette and scaling logic
  colors <- RColorBrewer::brewer.pal(11, "RdYlBu")
  heatmap_plot <- ggplot(long_matrix, aes(x = as.factor(Individual), y = SNP, fill = Importance)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("lightgrey", colors[6:11]),
      values = c(0, seq()),
      na.value = "lightgrey",
      limits = c(0, 1)
      # limits = c(0, max(non_zero_importance, na.rm = TRUE))
    ) +
    labs(
      title = paste(category, "Heterogeneity SNPs"),
      x = paste("Individuals (Ordered by ascending", ancestry_name, "%)"),
      y = "Genetic Variants (SNP)",
      fill = if (grepl("W", category)) "scaled positive\nW statistics" else "scaled selected\nW statistics"
    ) +
    scale_x_discrete(
      breaks = c("1", "5000", "10000", "15000", "20000"),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold", margin = margin(r = 8)),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
      legend.title = element_text(size = 8), 
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Save plot
  file_name <- file.path(save_path, paste0("fig3_", category, "_heatmap_sorted.png"))
  ggsave(file_name, plot = heatmap_plot, width = w, height = h, units = "in")
  print(paste(category, "heatmap saved to", file_name))
  return(heatmap_plot)
}



# main section: 
setwd("~/Project/Knockoff/run_real_data_julie")
snp_filepath <- "~/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv" # rows: samples, cols: PCs, EUR, ..., SNPs, Outcome(AD)
snp_data <- read.csv(snp_filepath)

eur <- snp_data$EUR
afr <- snp_data$AFR
amr <- snp_data$AMR
y <- snp_data$AD
non_eur <- 1 - eur

# TO CHANGE
exclRare <- TRUE # TO CHANGE
FDR <- 0.2 # TO CHANGE
ctrl_corr <- FALSE # TO CHANGE!!!!
load_directory = glue("~/Project/Knockoff/run_real_data_julie/run_standardize_TRUE_exclRare_TRUE_ctrlCorr_FALSE_0.9_fdr_0.2_model_Lasso_scale_pcs_update_w_compute")
genetic_variants <- snp_data[, grepl("chr", colnames(snp_data))] 
genetic_variants <- genetic_variants[, apply(genetic_variants, 2, var) != 0]
if (exclRare) {
  genetic_variants <- genetic_variants[, apply(genetic_variants, 2, sum) > dim(genetic_variants)[1]*2*0.05]
  genetic_variants <- genetic_variants[, apply(genetic_variants, 2, sum) < dim(genetic_variants)[1]*2*0.95]
}
if (ctrl_corr) {
  corr_mat <- cor(genetic_variants)
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
  genetic_variants <- genetic_variants[,representative_snps]
}
matrices <- load_matrices(load_directory)
colnames(matrices$scaled_selection_matrix) <- colnames(genetic_variants)
colnames(matrices$selection_matrix) <- colnames(genetic_variants)
colnames(matrices$W_statistic_matrix) <- colnames(genetic_variants)

selected_snp_names <- colnames(genetic_variants)[apply(matrices$selection_matrix, 2, any)]
print("Selected features based on the knockoff filter:")
print(selected_snp_names)
file_connection <- file(glue("{load_directory}/selected_snps.txt"), "w") # Open file for writing
cat(glue("Selected {length(selected_snp_names)} features based on the knockoff filter:"), file = file_connection, sep = "\n")
cat(selected_snp_names, file = file_connection, sep = "\n")
close(file_connection) # Close the file connection

######### plot heterogeneity based on scaled selection matrix ######### 
# select snps if it's selected in at least 10% of population
# selected_snp_names <- colnames(genetic_variants)[colMeans(matrices$selection_matrix) >= 0.5]

# Filter scaled_selection_matrix to include only selected SNPs
scaled_selection_matrix_filtered <- matrices$scaled_selection_matrix[, selected_snp_names, drop = FALSE]
max_non_zero <- max(scaled_selection_matrix_filtered)
min_non_zero <- min(scaled_selection_matrix_filtered[scaled_selection_matrix_filtered!=0])
# Categorize SNPs by non-zero counts using the filtered matrix
snp_categories <- categorize_snps_by_std_dev(scaled_selection_matrix_filtered)
print_snp_std_dev(snp_categories, scaled_selection_matrix_filtered)

low_snps <- select_top_snps(snp_categories$low_heterogeneity, "low", num_snps = 10)
medium_snps <- select_top_snps(snp_categories$medium_heterogeneity, "medium", num_snps = 10)
high_snps <- select_top_snps(snp_categories$high_heterogeneity, "high", num_snps = 10)

# low_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = low_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Low",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )
# 
# medium_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = medium_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Medium",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )
# 
# high_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = high_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "High",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )

# combined_W <- (low_heterogeneity_plot + theme(axis.text.y = element_blank())) + 
#   (medium_heterogeneity_plot + theme(axis.text.y = element_blank())) + 
#   (high_heterogeneity_plot + theme(axis.text.y = element_blank())) +
#   plot_layout(guides = "collect")
# ggsave(glue("{load_directory}/for_poster.png"), combined_W, width = 10, height = 2.8, dpi = 400)
# all_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = matrices$scaled_selection_matrix,
#   selected_snps = selected_snp_names,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "ALL",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 8,
#   h = 5
# )

# scaled
all_heterogeneity_plot <- plot_heatmap(
  scaled_selection_matrix = matrices$scaled_selection_matrix,
  selected_snps = snp_categories$ranked_snps, # Here give ranked snps
  ancestry = non_eur,
  ancestry_name = "non-EUR",
  category = "ranked ALL",
  save_path = load_directory,
  max_non_zero = max_non_zero,
  min_non_zero = min_non_zero,
  w = 8,
  h = 5
)
# combined_W <- high_heterogeneity_plot / medium_heterogeneity_plot / low_heterogeneity_plot
# 
# ggsave(glue("{load_directory}/combined_scaled_selected_lmh.png"), combined_W, width = 8, height = 6, units = "in", dpi = 400)

#################################################
######### plot heterogeneity based on W ######### 
#################################################
W = matrices$W_statistic_matrix

# Replace 0s with NA to avoid division by zero
max_vals <- apply(W, 2, max, na.rm = TRUE)
safe_max <- ifelse(max_vals == 0, NA, max_vals)
scaled_W <- sweep(W, 2, safe_max, "/")
# Replace columns with NA (i.e., originally 0 max) with original W
scaled_W[, is.na(safe_max)] <- W[, is.na(safe_max)]
scaled_W[scaled_W < 0] <- 0


# scaled_W_2 <- t((t(W) / apply(W, 2, max, na.rm = TRUE)))
scaled_W_selected <- scaled_W[, selected_snp_names, drop = FALSE]
max_non_zero <- max(scaled_W_selected)
min_non_zero <- min(scaled_W_selected[scaled_W_selected>0])

snp_categories <- categorize_snps_by_std_dev(scaled_W_selected) # var calculated include negative W
print_snp_std_dev(snp_categories, scaled_W_selected)

low_snps <- select_top_snps(snp_categories$low_heterogeneity, "low", num_snps = 10)
medium_snps <- select_top_snps(snp_categories$medium_heterogeneity, "medium", num_snps = 10)
high_snps <- select_top_snps(snp_categories$high_heterogeneity, "high", num_snps = 10)

# low_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = scaled_W,
#   selected_snps = low_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Low-W",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )
# 
# medium_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = scaled_W,
#   selected_snps = medium_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "Medium-W",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )
# 
# high_heterogeneity_plot <- plot_heatmap(
#   scaled_selection_matrix = scaled_W,
#   selected_snps = high_snps,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "High-W",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 7,
#   h = 2
# )
# 
# all_W_heterogeneity_plot_scaled <- plot_heatmap(
#   scaled_selection_matrix = scaled_W,
#   selected_snps = selected_snp_names,
#   ancestry = non_eur,
#   ancestry_name = "non-EUR",
#   category = "W_all_scaled",
#   save_path = load_directory,
#   max_non_zero = max_non_zero,
#   min_non_zero = min_non_zero,
#   w = 8,
#   h = 5
# )

# scaled
all_ranked_W_heterogeneity_plot_scaled <- plot_heatmap(
  scaled_selection_matrix = scaled_W,
  selected_snps = snp_categories$ranked_snps, # Here give ranked snps
  ancestry = non_eur,
  ancestry_name = "non-EUR",
  category = "ranked_W_all_scaled",
  save_path = load_directory,
  max_non_zero = max_non_zero,
  min_non_zero = min_non_zero,
  w = 8,
  h = 5
)

# combined low medium high
# combined_W <- high_heterogeneity_plot / medium_heterogeneity_plot / low_heterogeneity_plot
#             
# ggsave(glue("{load_directory}/combined_W_lmh.png"), combined_W, width = 8, height = 6, units = "in", dpi = 400)

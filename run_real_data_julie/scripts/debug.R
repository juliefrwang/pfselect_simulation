

knockoff_filter <- function(local_feature_importance, local_feature_importance_knockoff, FDR_rate = 0.1) {
  W <- local_feature_importance - local_feature_importance_knockoff
  q_values <- matrix(NA, nrow = nrow(W), ncol = ncol(W))  # Each element represents a q-value for individual i and SNP j
  
  # Calculate q-values for each individual across all SNPs
  for (i in 1:nrow(W)) {
    # Calculate MK statistics for individual i across all SNPs
    MK_stat <- MK.statistic(local_feature_importance[i, ], local_feature_importance_knockoff[i, ])
    q <- MK.q.byStat(MK_stat[,'kappa'], MK_stat[,'tau'], M=1, Rej.Bound = 10000)
    q_values[i, ] <- q
  }
  # Create the selection matrix S_ij: S_ij = I(q_ij <= FDR_rate)
  S_ij <- q_values <= FDR_rate
  
  # Create scaled selection matrix using W matrix and selection matrix S_ij
  scaled_selection_matrix <- t((t(W) / apply(W, 2, max, na.rm = TRUE))) * S_ij
  
  return(list(S_ij = S_ij, scaled_selection_matrix = scaled_selection_matrix, W = W, q_values = q_values))
}

W1 <- read.csv("~/Project/Knockoff/Lab_RY/lab_1_6/results/W_statistic_matrix.csv", row.names = 1)
W2 <- read.csv("~/Project/Knockoff/run_real_data_julie/runs/W_statistic_matrix.csv",)
W3 <- read.csv("~/Project/Knockoff/run_real_data_julie/runs_standardize_true_two_pcs/W_statistic_matrix.csv",)

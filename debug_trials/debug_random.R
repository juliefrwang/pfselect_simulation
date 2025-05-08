
set.seed(130)
snps <- matrix(rnorm(5000 * 20), nrow = 5000, ncol = 20)


# ancestry
# ancestry <- as.matrix(runif(5000))
ancestry <- as.matrix(c(rep(0, 5000)))

# target variable
true_beta <- numeric(20)
true_beta[1:5] <- c(0.1, 0.1, -0.1, 0.1, -0.1)

true_y <- numeric(nrow(snps))

for (i in 1:nrow(snps)) {
  beta_mod <- true_beta
  #beta_mod[1] <- beta_mod[1] * ancestry[i]
  #beta_mod[2] <- beta_mod[2] * ancestry[i]
  true_y[i] <- snps[i,] %*% beta_mod
}

#hist(true_y)

true_y <- true_y + rnorm(nrow(snps))


# Initialize result storage
num_rows <- nrow(snps)  # Number of individuals
fdp_list <- numeric(num_rows)
power_list <- numeric(num_rows)

# Define the indices of significant beta valuesr
sig_beta_inx <- which(true_beta!= 0)


# start replicate
num_replicate = 3
for (i in 1:num_replicate) {
  snps_knockoff <- generate_knockoff_data(snps)
  print(snps_knockoff[1,1:5])
  results <- get_importance_matrices(
    genetic_variants = snps,
    genetic_variants_knockoff = snps_knockoff,
    additional_covariates = NULL,
    Z = ancestry, # or other variables like 'pcs'
    y = true_y,
    n_folds = 5,
    FDR_rate = 0.1
  )
  
  # Loop over all individuals
  for (i in 1:num_rows) {
    selected <- which(results$selection_matrix[i, ] == 1)
    
    # Calculate true and false positives
    num_true_positive <- sum(selected %in% sig_beta_inx)
    num_false_positive <- length(selected) - num_true_positive
    
    # Compute FDP and power
    FDP <- num_false_positive / max(1, length(selected))
    power <- num_true_positive / max(1, length(sig_beta_inx))
    
    # Store results
    fdp_list[i] <- fdp_list[i] + FDP
    power_list[i] <- power_list[i] + power
  }
}

fdp_list_avg <- fdp_list / num_replicate
power_list_avg <- power_list / num_replicate

ordered_indices <- order(ancestry)  # Ascending order of v
fdp_list_ordered <- fdp_list_avg[ordered_indices]
power_list_ordered <- power_list_avg[ordered_indices]

plot(c(1:num_rows), fdp_list_ordered, type = "b", col = "red", pch = 19, lwd = 2, ylim = c(0, 1.05))
lines(c(1:num_rows), power_list_ordered, type = "b", col = "blue", pch = 19, lwd = 2)
abline(h = 0.1, col = "green", lwd = 2, lty = 2)


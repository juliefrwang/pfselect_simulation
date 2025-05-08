library(adelie)
library(glmnet)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(patchwork)
library(gglasso)
library(glue)
# sanity check of lasso and group lasso 
generate_knockoff_data <- function(genetic_variants) {
  # Extract only the SNP columns (those that start with "chr") to create the genotype matrix
  G <- as.matrix(genetic_variants)
  
  # Basic operations for knockoff preparation
  sd.G <- apply(G, 2, sd)
  cor.G <- matrix(as.numeric(corpcor::cor.shrink(G, verbose = FALSE)), nrow = ncol(G))
  
  # Generate preliminary fit for knockoff creation
  fit.prelim.M.sdp <- GhostKnockoff.prelim(cor.G, M = 1, method = 'sdp')
  
  
  # Second order knockoffs - sdp
  mu.G <- colMeans(G)
  scale.G <- t((t(G) - mu.G) / sd.G)  # Center and scale each SNP column
  # Generate the knockoff matrix
  M <- 1
  E <- matrix(rnorm(nrow(scale.G) * M * ncol(scale.G)), nrow(scale.G), M * ncol(scale.G))
  # Convert V.left to a standard matrix as it’s matrix-like from the Matrix package
  V.left <- as.matrix(fit.prelim.M.sdp$V.left)
  scale.G_k1.M <- t(apply(scale.G %*% t(fit.prelim.M.sdp$P.each), 1, rep, times = M)) + E %*% t(V.left)
  
  # Transform back to original scale
  temp.G_k1.M <- t(rep(sd.G, M) * t(scale.G_k1.M) + rep(mu.G, M))
  G_k <- as.matrix(temp.G_k1.M)  # Knockoff matrix
  
  print(paste("Genotype knockoff dimensions:", dim(G_k)[1], "samples x", dim(G_k)[2], "SNPs"))
  
  return(G_k)
}

set_penalty_factors <- function(num_snps, num_z, num_u) {
  penalty_factors <- c(rep(0, num_u), rep(1, 2 * num_snps + 2 * num_snps * num_z))
  return(penalty_factors)
}

set_penalty_factors_grpLasso <- function(num_snps, num_z, num_u) {
  # Default is square-root of group sizes for each group.
  penalty_for_grp <- sqrt(1+num_z) # penalty for groups except unpenalized covariates
  num_grp <- 2 * num_snps # number of groups except unpenalized covariates
  penalty_factors <- c(rep(0, num_u), rep(penalty_for_grp, num_grp))
  return(penalty_factors)
}


extract_coefficients <- function(cv_fit) {
  coefs <- coef(cv_fit, s = "lambda.min")
  print(paste("Number of coefficients:", length(coefs)))
  return(coefs)
}

extract_grp_lasso_coefficients <- function(cv_fit, recovery_index) {
  
  # Note that if the default standardize=TRUE was used in fitting the grpnet 
  # object, the coefficients reported are for the standardized inputs. However, 
  # the predict function will apply the stored standardization to newx and give 
  # the correct predictions.
  
  ordered_coefs <- coef(cv_fit, lambda = "lambda.1se")$betas 
  coefs <- ordered_coefs[recovery_index]
  print(paste("Number of coefficients:", length(coefs)))
  return(coefs)
}

# simulated data
sim_path <- "seed_379_b_0.35_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_FALSE_standardize_FALSE_w_AMR_sampleSize_not_equal"
simulated_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
pcs <- simulated_data[,c("PC1", "PC2", "PC3", "PC4")]
snps <- simulated_data[,grepl("chr", colnames(simulated_data))]
snps_knockoff <- generate_knockoff_data(snps)
# X <- as.matrix(cbind(snps, snps_knockoff))

######################### 1. simple case. ######################### 
set.seed(100)
# snps <- matrix(rnorm(5000*124), nrow = 5000, ncol = 124)
non_zero_beta <- rnorm(ncol(snps) %/% 2)
zeros <- rep(0, ncol(snps) %/% 2)
beta <- as.vector(c(non_zero_beta, zeros))
if (length(beta) != ncol(snps)) {
  beta <- c(beta,0)
}
y <- as.matrix(snps) %*% beta + 2 + rnorm(nrow(snps))

lasso_group <- seq(1,ncol(snps))
glasso_group <- rep(1:(ncol(snps)/2), each=2)
if (length(glasso_group) != ncol(snps)) {
  glasso_group <- c(beta,glasso_group[length(glasso_group)])
}

glmnet_lasso_cv_fit <- cv.glmnet(x = as.matrix(snps), 
                    y = y, 
                    family = "gaussian", 
                    alpha = 1, 
                    n_folds = 10, 
                    standardize = TRUE)

grpnet_lasso_cv_fit <- cv.grpnet(X = as.matrix(snps), 
                    glm = glm.gaussian(y), 
                    alpha = 1, 
                    n_folds = 10, 
                    adev_tol = 0.95,
                    ddev_tol = 1e-5,
                    lmda_path_size = 150,
                    standardize = TRUE)


gglasso_lasso_cv_fit <- cv.gglasso(x = as.matrix(snps), 
                                   y = y, 
                                   group = lasso_group,
                                   loss = "ls",
                                   pred.loss = "L2", # mean square error used by least squares
                                   nfolds = 10)

gglasso_lasso_cv_fit_3 <- cv.gglasso(x = as.matrix(snps), 
                                     y = y, 
                                     group = glasso_group,
                                     loss = "ls",
                                     pred.loss = "L2", # mean square error used by least squares
                                     nfolds = 10)
# print(glmnet_lasso_cv_fit)
# print(grpnet_lasso_cv_fit)
# print(gglasso_lasso_cv_fit)

plot(glmnet_lasso_cv_fit)
plot(grpnet_lasso_cv_fit)
plot(gglasso_lasso_cv_fit)


# Mannually check the CV path and compute the MSE
library(adelie)
library(caret)  # for folds

set.seed(123)
folds <- createFolds(y, k = 10)
lambda_seq_1 <- glmnet_lasso_cv_fit$lambda
lambda_seq_2 <- grpnet_lasso_cv_fit$lambda

fold_mse_1 <- matrix(NA, nrow = length(lambda_seq_1), ncol = 10)
fold_mse_2 <- matrix(NA, nrow = length(lambda_seq_2), ncol = 10)
x = as.matrix(snps)

for (i in 1:10) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(y), test_idx)
  
  x_train <- x[train_idx, ]
  y_train <- y[train_idx]
  x_test  <- x[test_idx, ]
  y_test  <- y[test_idx]
  
  # # glmnet fit
  # glmnet.control(devmax = 0.99, fdev = 1e-5)
  # glmnet_fit <- glmnet(x_train, y_train, alpha = 1, family = "gaussian", 
  #                      lambda = lambda_seq_1, standardize = TRUE)
  # preds1 <- predict(glmnet_fit, newx = x_test, lambda = lambda_seq_1, type = "response") # matrix: n_test x length(lambda_seq)
  # 
  # # Compute MSE for each lambda
  # for (j in seq_along(lambda_seq_1)) {
  #   y_pred <- preds1[, j]
  #   fold_mse_1[j, i] <- mean((y_test - y_pred)^2)
  # }
  # 
  # grpnet fit
  grpnet_fit <- grpnet(x_train, glm.gaussian(y_train), alpha = 1, lambda = lambda_seq_2, adev_tol = 0.95, ddev_tol = 1e-5, standardize = TRUE)
  
  preds2 <- predict(grpnet_fit, newx = x_test, lambda = lambda_seq_2, type = "response")  # matrix: n_test x length(lambda_seq)
  
  # Compute MSE for each lambda
  for (j in seq_along(lambda_seq_2)) {
    y_pred <- preds2[, j]
    fold_mse_2[j, i] <- mean((y_test - y_pred)^2)
  }
}

mean_mse_per_lambda_glmnet <- rowMeans(fold_mse_1)
# Average MSE across folds for each lambda
mean_mse_per_lambda_grpnet <- rowMeans(fold_mse_2)





#### fit Lasso with regularization path ####
set.seed(100)
# snps <- matrix(rnorm(5000*124), nrow = 5000, ncol = 124)
non_zero_beta <- rnorm(ncol(snps) %/% 2)
zeros <- rep(0, ncol(snps) %/% 2)
beta <- as.vector(c(non_zero_beta, zeros))
if (length(beta) != ncol(snps)) {
  beta <- c(beta,0)
}
y <- as.matrix(snps) %*% beta + 2 + rnorm(nrow(snps))
glasso_group <- seq(1,ncol(snps))

################ 1. glmnet ################   
# ?glmnet.control
# glmnet.measures()
glmnet.control(devmax = 0.95, fdev = 1e-5)
glmnet_fit <- glmnet(snps, y, alpha = 1, nlambda = 100, standardize = TRUE)
plot(glmnet_fit, xvar = "lambda")
df1 <- print(glmnet_fit)
# plot(glmnet_fit, xvar = "dev")

################ 2. grpnet ################  
grpnet_fit <- grpnet(snps, glm.gaussian(y), alpha = 1, lambda =lambda_seq_1, adev_tol = 0.95, ddev_tol = 1e-5, lmda_path_size = 100, standardize = TRUE)
df2 <- print(grpnet_fit)
plot(grpnet_fit)

glmnet_fit <- glmnet(snps, y, alpha = 1, lambda=0.2, standardize = TRUE)
glmnet_fit$beta

grpnet_fit <- grpnet(snps, glm.gaussian(y), alpha = 1, lambda =0.2, standardize = TRUE)


gglasso_fit <- gglasso(as.matrix(snps), y, group = lasso_group, lambda=0.2)
gglasso_fit$beta

a<-cbind(glmnet_fit$beta,t(coef(grpnet_fit)$betas),gglasso_fit$beta)



################ 3. gglasso ################  
# https://cran.r-project.org/web/packages/gglasso/gglasso.pdf
gglasso_fit <- gglasso(as.matrix(snps), y, group = lasso_group, loss = "ls")
df3 <- print(gglasso_fit)
plot(gglasso_fit)

# cbind together and show
pad_df <- function(df, n) {
  if (nrow(df) < n) {
    padding <- as.data.frame(matrix(NA, nrow = n - nrow(df), ncol = ncol(df)))
    colnames(padding) <- colnames(df)
    df <- rbind(df, padding)
  }
  return(df)
}


df1 <- pad_df(df1, 100)
df4 <- cbind(df1, df2, df3)

coefs <- coef(lasso_cv_fit, s = "lambda.1se")
coefs <- coef(glasso_cv_fit, s = "lambda.1se")
# fix lambda
glmnet_fit_2 <- glmnet(snps, y, alpha = 1, lambda = seq(0.1,0.01,-0.001), nlambda = 100, standardize = TRUE)
grpnet_fit_2 <- grpnet(snps, glm.gaussian(y), alpha = 1, lambda = seq(0.1,0.01,-0.001), adev_tol = 0.95, ddev_tol = 0, lmda_path_size = 100, standardize = TRUE)
df1.2 <- print(glmnet_fit_2)
df2.2 <- print(grpnet_fit_2)
df3.2 <- cbind(df1.2, df2.2)




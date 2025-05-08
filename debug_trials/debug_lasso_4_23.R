library(adelie)
library(glmnet)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(patchwork)
library(gglasso)
library(glue)


set.seed(100)

# random x
# x <- matrix(rnorm(5000*124), nrow = 5000, ncol = 124)

# real data
sim_path <- "seed_379_b_0.35_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_FALSE_standardize_FALSE_w_AMR_sampleSize_not_equal"
simulated_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
snps <- simulated_data[,grepl("chr", colnames(simulated_data))]
# ctrl_corr = TRUE
# if (ctrl_corr) {
#   corr_mat <- cor(snps)
#   Sigma.distance = as.dist(1 - abs(corr_mat))
#   fit = hclust(Sigma.distance, method="single")
#   corr_max <- 0.9
#   clusters = cutree(fit, h=1-corr_max)
#   representative_snps <- numeric(length(unique(clusters)))
# 
#   # Find representative SNP for each cluster
#   for (k in unique(clusters)) {
#     cluster_indices <- which(clusters == k)
#     if (length(cluster_indices) == 1) {
#       representative_snps[k] <- cluster_indices[1]
#     } else {
#       # Calculate sum of absolute correlations within the cluster
#       cluster_corr_sum <- colSums(abs(corr_mat)[cluster_indices, cluster_indices])
#       # Find the SNP with the highest sum of correlations
#       representative_snps[k] <- cluster_indices[which.max(cluster_corr_sum)]
#     }
#   }
#   snps <- snps[,representative_snps]
# }
x = as.matrix(snps)
non_zero_beta <- rnorm(ncol(x) %/% 2)
zeros <- rep(0, ncol(x) %/% 2)
beta <- as.vector(c(non_zero_beta, zeros))
if (length(beta) != ncol(x)) {
  beta <- c(beta,0)
}
y <- 2 + x %*% beta + rnorm(nrow(x))

# for gglasso input
lasso_group <- seq(1,ncol(x))

# standardize
x <- matrix(rnorm(3*2), nrow = 3, ncol = 2)
x_mean<-colMeans(x)
x_sd<-apply(x,2,sd)
x_standard<-t((t(x)-x_mean)/x_sd)

# ################ 1. glmnet ################   
glmnet_fit <- glmnet(x_standard, y, alpha = 1, 
                     nlambda = 100, 
                     lambda.min.ratio = 0.0001, 
                     standardize = TRUE)
plot(glmnet_fit, xvar = "lambda")
log(glmnet_fit$lambda)[2:10]-log(glmnet_fit$lambda)[1:9]
################ 2. grpnet ################
grpnet_fit <- grpnet(x_standard, glm.gaussian(y), alpha = 1, 
                     lmda_path_size = 100, 
                     min_ratio = 0.0001,
                     standardize = TRUE,
                     adev_tol = 0.999,
                     ddev_tol = 1e-5)

plot(grpnet_fit)
ll<-print(grpnet_fit)
log(ll$Lambda)[2:10]-log(ll$Lambda)[1:9]

################ 3. gglasso ################
# https://cran.r-project.org/web/packages/gglasso/gglasso.pdf
gglasso_fit <- gglasso(x_standard, y, group = lasso_group, lambda.factor=0.0001)
plot(gglasso_fit)
log(gglasso_fit$lambda)[2:10]-log(gglasso_fit$lambda)[1:9]

################ fix lambda #################
glmnet_fix_lambda <- glmnet(x, y, alpha = 1, lambda=0.2, standardize = TRUE)
grpnet_fix_lambda <- grpnet(x, glm.gaussian(y), alpha = 1, lambda =0.2, 
                            min_ratio = 0.0001, 
                            standardize = TRUE,
                            adev_tol = 0.999,
                            ddev_tol = 1e-5)
gglasso_fix_lambda <- gglasso(x, y, group = lasso_group, lambda=0.2,lambda.factor=0.0001)

df <-cbind(glmnet_fix_lambda$beta, t(coef(grpnet_fix_lambda)$betas), gglasso_fix_lambda$beta)
#############
glmnet_fix_lambda <- glmnet(x_standard, y, alpha = 1, lambda=0.2, standardize = TRUE)
grpnet_fix_lambda <- grpnet(x_standard, glm.gaussian(y), alpha = 1, lambda =0.2, 
                            min_ratio = 0.0001, 
                            standardize = TRUE,
                            adev_tol = 0.999,
                            ddev_tol = 1e-5)
gglasso_fix_lambda <- gglasso(x_standard, y, group = lasso_group, lambda=0.2,lambda.factor=0.0001)

df_standard <-cbind(glmnet_fix_lambda$beta, t(coef(grpnet_fix_lambda)$betas), gglasso_fix_lambda$beta)


df_all<-as.matrix(round(cbind(df,df_standard),4))


############## cross validation ###########
glmnet_lasso_cv_fit <- cv.glmnet(x = x_standard, 
                                 y = y, 
                                 family = "gaussian", 
                                 alpha = 1, 
                                 n_folds = 10, 
                                 standardize = TRUE)

grpnet_lasso_cv_fit <- cv.grpnet(X = x_standard, 
                                 glm = glm.gaussian(y), 
                                 alpha = 1, 
                                 n_folds = 10, 
                                 adev_tol = 0.999,
                                 ddev_tol = 1e-5,
                                 lmda_path_size = 100,
                                 min_ratio = 0.0001,
                                 standardize = TRUE)

gglasso_lasso_cv_fit <- cv.gglasso(x = x_standard, 
                                   y = y, 
                                   group = lasso_group,
                                   loss = "ls",
                                   lambda.factor=0.0001,
                                   pred.loss = "L2", # mean square error used by least squares
                                   nfolds = 10)

print(glmnet_lasso_cv_fit)
print(grpnet_lasso_cv_fit)
print(gglasso_lasso_cv_fit)

plot(glmnet_lasso_cv_fit)
plot(grpnet_lasso_cv_fit)
plot(gglasso_lasso_cv_fit)


# standardize = FALSE
# ################ 1. glmnet ################   
glmnet_fit_no_scale <- glmnet(x, y, alpha = 1, 
                     nlambda = 100, 
                     lambda.min.ratio = 0.0001, 
                     standardize = FALSE)
plot(glmnet_fit_no_scale)

log(glmnet_fit_no_scale$lambda)[2:10]-log(glmnet_fit_no_scale$lambda)[1:9]

################ 2. grpnet ################
grpnet_fit_no_scale <- grpnet(x, glm.gaussian(y), alpha = 1, 
                     lmda_path_size = 100, 
                     min_ratio = 0.0001,
                     standardize = FALSE,
                     adev_tol = 0.999,
                     ddev_tol = 1e-5)
plot(grpnet_fit_no_scale)
ll_no_scale<-print(grpnet_fit_no_scale)
log(ll_no_scale$Lambda)[2:10]-log(ll_no_scale$Lambda)[1:9]
##### cross validate with standardize = FALSE 
glmnet_lasso_cv_fit_no_scale <- cv.glmnet(x = x, 
                                 y = y, 
                                 family = "gaussian", 
                                 alpha = 1, 
                                 n_folds = 10, 
                                 standardize = FALSE)

grpnet_lasso_cv_fit_no_scale <- cv.grpnet(X = x, 
                                 glm = glm.gaussian(y), 
                                 alpha = 1, 
                                 n_folds = 10, 
                                 adev_tol = 0.999,
                                 ddev_tol = 1e-5,
                                 lmda_path_size = 100,
                                 min_ratio = 0.0001,
                                 standardize = FALSE)

library(adelie)
library(glmnet)
library(jsonlite)
library(ggplot2)
library(reshape2)
library(patchwork)

# check LASSO whether 
sim_path <- "seed_379_b_0.35_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_sampleSize_not_equal_model_Lasso_no_amr"
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
simulated_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))
snps <- simulated_data[, grepl("chr", colnames(simulated_data))]
pcs = simulated_data[,c("PC1","PC2","PC3","PC4")]
beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
res1.lasso <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/rep_70.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx
amr_idx <- beta$amr_idx
num_snps <- ncol(snps)

# check coef
for(x in afr_idx) {
  ind_all_x = 1 + 4 + seq(x, x + 9 * num_snps, by = num_snps)
  print(unlist(res1.lasso$coefs[ind_all_x]))
}

# check corra
corr_mat = as.matrix(cor(snps))
corr_mat[,afr_idx[6]][corr_mat[,afr_idx[6]] > 0.75]
corr_mat[,eur_idx[1]][corr_mat[,eur_idx[1]] > 0.75]


# Check W statistics
head(res1.lasso$W_statistic_matrix[,eur_idx])
# compute W
x <- eur_idx[8] # which snp
ind_all_x = 1 + 4 + seq(x, x + 9 * num_snps, by = num_snps)
coefs_lasso = unlist(res1.lasso$coefs[ind_all_x])
i = 10 # which individual
pcs = simulated_data[i, c("PC1","PC2","PC3","PC4")]
T_origin = abs(coefs_lasso[1] + sum(coefs_lasso[2:5] * pcs)) 
T_tilda = abs(coefs_lasso[6] + sum(coefs_lasso[7:10] * pcs))
print(glue("coefs: {coefs_lasso}"))
print(glue("T_origin: {T_origin}"))
print(glue("T_tilda: {T_tilda}"))
print(T_origin - T_tilda)
print(glue("W: {res1.lasso$W_statistic_matrix[i,x]}"))
print(list(non_zero = sum(res1.lasso$coefs != 0), 
           total = length(res1.lasso$coefs)))
print(glue("predicted y: {sum(simulated_data[i, 1:124] * res1.lasso$coefs)}"))
print(glue("true y: {simulated_data[i,]$true_y_no_noise}"))

# check Group LASSO whether 
sim_path <- ("seed_379_b_0.35_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_sampleSize_not_equal_model_Grplasso_no_amr_1se")
simulation_results <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulation_result.csv"))
simulated_data <- read.csv(glue("~/Project/Knockoff/simulation_res/{sim_path}/simulated_true_data.csv"))

beta <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/beta_idx.json"))
res1.grplasso <- fromJSON(glue("~/Project/Knockoff/simulation_res/{sim_path}/Wstats/rep_11.json"))
true_idx <- beta$true_idx
homo_idx <- beta$homo_idx
eur_idx <- beta$eur_idx
afr_idx <- beta$afr_idx
amr_idx <- beta$amr_idx
num_snps <- beta$num_snps

for(x in eur_idx) {
  ind_all_x = 4 + seq(x, x + 9 * num_snps, by = num_snps)
  print(unlist(res1.grplasso$coefs[ind_all_x]))
}

for(x in afr_idx) {
  ind_all_x = 1 + 4 + seq(x, x + 9 * num_snps, by = num_snps)
  print(unlist(res1.grplasso$coefs[ind_all_x]))
}
print(res1.grplasso$coefs[1:4])

# Check W statistics
head(res1.grplasso$W_statistic_matrix[,eur_idx])
# compute W
x <- afr_idx[6] # which snp
ind_all_x = 4 + seq(x, x + 9 * num_snps, by = num_snps)
coefs_grplasso = unlist(res1.grplasso$coefs[ind_all_x])
i = 10 # which individual
pcs = simulated_data[i, c("PC1","PC2","PC3","PC4")]
T_origin = abs(coefs_grplasso[1] + sum(coefs_grplasso[2:5] * pcs)) 
T_tilda = abs(coefs_grplasso[6] + sum(coefs_grplasso[7:10] * pcs))
print(glue("coefs: {coefs_grplasso}"))
print(glue("T_origin: {T_origin}"))
print(glue("T_tilda: {T_tilda}"))
print(T_origin - T_tilda)
print(glue("W: {res1.grplasso$W_statistic_matrix[i,x]}"))
print(list(non_zero = sum(res1.grplasso$coefs != 0), 
           total = length(res1.grplasso$coefs)))

print(glue("predicted y: {sum(simulated_data[i, 1:124] * res1$coefs)}"))
print(glue("true y: {simulated_data[i,]$true_y_no_noise}"))

# Side-by-Side Barplot
coefs_lasso = unlist(res1.lasso$coefs)[-1,1]
coefs_grplasso = unlist(res1.grplasso$coefs)[,1]
df <- data.frame(
  Index = 1:length(coefs_lasso),
  Lasso = coefs_lasso,
  GroupLasso = coefs_grplasso
)
# Filter only non-zero in either
df_nonzero <- subset(df, Lasso != 0 | GroupLasso != 0)

# Plot both on the same axes
lasso_plt <- ggplot(subset(df, Lasso != 0), aes(x = Index, y = Lasso, color="Lasso")) +
  geom_point(shape = 16, size = 2) +
  scale_color_manual(values = c("Lasso" = "blue", "GroupLasso" = "tomato")) +
  theme_minimal() +
  labs(title = "non-zero Lasso Coefficients", x = "Index", y = "Coefficient Value", color = "Model")+
  ylim(c(-26, 26)) + 
  xlim(c(0,1245))

grplasso_plt <- ggplot(subset(df, GroupLasso != 0), aes(x = Index, y = GroupLasso, color="GroupLasso")) +
  geom_point(shape = 16, size = 2) +
  scale_color_manual(values = c("Lasso" = "blue", "GroupLasso" = "tomato")) +
  theme_minimal() +
  labs(title = "non-zero group Lasso Coefficients", x = "Index", y = "Coefficient Value", color = "Model")+
  ylim(c(-26, 26)) + 
  xlim(c(0,1245))

lasso_plt + grplasso_plt + plot_layout(guides = "collect") + plot_annotation(title = "Non-Zero Coefficients: Lasso vs Group Lasso")

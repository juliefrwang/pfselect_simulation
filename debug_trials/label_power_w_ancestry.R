library(ggplot2)
library(patchwork)
library(glue)
library(jsonlite)
library(parallel)
library(dplyr)

# get pcs
set.seed(379)
snp_filepath <- "/Users/juliew/Project/Knockoff/Lab_RY/lab_1_6/data/final_snp_data.csv"
snp_data <- read.csv(snp_filepath)
snp_data <- snp_data %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 

# Generate a 5000 sample
sample_data <- snp_data[sample(1:nrow(snp_data), 5000),] 
# pcs and ancestry
ancestry <- as.matrix(sample_data[,c("SAS","EAS","AMR","AFR","EUR")])
ans <- sample_data[,"ans_label"]

# read output
sim_path <- "seed_379_b_0.12_FDR_0.2_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_TRUE_sampleSize_not_equal_Lasso_ctrlAnsMaf_FALSE_scale_pcs_varied_b"
setwd(glue("~/Project/Knockoff/simulation_res/{sim_path}"))
df <- read.csv(glue("simulation_result.csv"))

fdp_list <- df$fdp_list
#power_list <- df$power_list
power_eur_list <- df$power_eur_list
power_afr_list <- df$power_afr_list
power_homo_list <- df$power_homo_list

# prepare plot
ordered_indices <- order(ancestry[,'EUR'])  # Ascending order of v
fdp_list_ordered <- fdp_list[ordered_indices]
power_eur_list_ordered <- power_eur_list[ordered_indices]
power_homo_list_ordered <- power_homo_list[ordered_indices]
ans_ordered <- ans[ordered_indices]

# create top panel: panel vs x (ordered by EUR)
set2_colors <- brewer.pal(n = 8, name = "Set2")
ancestry_colors <- c("SAS" = set2_colors[1],   
                     "EAS" = set2_colors[2],  
                     "AMR" = set2_colors[3],  
                     "AFR" = set2_colors[4],  
                     "EUR" = set2_colors[5],  
                     "Admixed" = set2_colors[8]) 

power_plot_1 <- ggplot() +
  geom_point(aes(x = c(1:num_rows), y = power_homo_list_ordered, color = ans_ordered), 
             alpha = 0.7, size = 1) + 
  # geom_point(aes(x = c(1:num_rows), y = power_homo_list_ordered), 
  #            alpha = 0.7, size = 1, color = 'orange') + 
  ylim(0,1.0) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = "Power", color = NULL) +
  scale_colour_manual(values = ancestry_colors) 

ggsave("labeled.png", power_plot_1, width = 6, height = 4)


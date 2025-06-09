df <- read.csv("~/Project/Knockoff/Lab_RY/lab_1_6/data/Ancestry_ADSP_WGS_36k_geno0.05_SNPWeight_Yann.tsv", sep = '\t')
pcs <- df[c("PC1", "PC2", "PC3", "PC4")]
df <- df %>%
  mutate(ans_label = ifelse(EUR> 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed"))))))
for_plot <- data.frame(
  PC1 <- df$PC1,
  PC2 <- df$PC2,
  PC3 <- df$PC3,
  PC4 <- df$PC4,
  ans_label <- df$ans_label
)
set2_colors <- brewer.pal(n = 8, name = "Set2")
ancestry_colors <- c("SAS" = set2_colors[1],   
                     "EAS" = set2_colors[2],  
                     "AMR" = set2_colors[3],  
                     "AFR" = set2_colors[4],  
                     "EUR" = set2_colors[5],  
                     "Admixed" = set2_colors[8]) 

# Map the colors to the ancestry labels

# pc1 and pc2
p1 <- ggplot(for_plot, aes(x = PC1, y = PC2, color = ans_label)) +
  geom_point(alpha = .7) +
  labs(title = "ADSP Cohort",
       x = "PC1",
       y = "PC2",
       color = "Ancestry label") +
  base_theme +
  scale_colour_manual(values = ancestry_colors) +
  coord_fixed(ratio = 1)

# pc3 and pc4
p2 <- ggplot(for_plot, aes(x = PC3, y = PC4, color = ans_label)) +
  geom_point(alpha = .7) +
  labs(title = "ADSP Cohort",
       x = "PC3",
       y = "PC4",
       color = "Ancestry label") +
  base_theme +
  scale_colour_manual(values = ancestry_colors) +
  coord_fixed(ratio = 1)

combined <- p1 + p2
combined

ggsave("~/Project/Knockoff/Lab_RY/lab_1_6/data/pcs.png", combined, width=8, height=5)

pc_scale <- scale(pcs)

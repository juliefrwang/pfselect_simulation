library(PFSelectcopy)
library(glue)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(patchwork)

which.path <- "seed_379_b_0.4_rep_100_true_40_hete_30_exclRare_TRUE_ctrlCorr_TRUE_standardize_FALSE_w_AMRbeta_exclAMRsample_sampleSize_equal"
setwd(glue("~/Project/Knockoff/simulation_res/{which.path}")) #TOCHANGE
df.true <- read.csv("simulated_true_data.csv")
df.res <- read.csv("simulation_result.csv")
fdp_list <- df.res$fdp_list
power_list <- df.res$power_list
power_homo_list <- df.res$power_homo_list

df.true <- df.true %>%
  mutate(ans_label = ifelse(EUR > 0.75, "EUR",
                     ifelse(AFR > 0.75, "AFR",
                     ifelse(AMR > 0.75, "AMR",
                     ifelse(EAS > 0.75, "EAS",
                     ifelse(SAS > 0.75, "SAS", "Admixed")))))) 
# plot overall beta along different ancestry label
for.plot <- data.frame(power = power_list, ans_label = df.true$ans_label)
for.plot$ans_label <- factor(for.plot$ans_label, levels = c("Admixed", "AFR", "AMR", "EUR"))
for.plot <- for.plot %>%
  group_by(ans_label) %>%
  mutate(original_order = row_number()) %>%
  ungroup() %>%
  arrange(ans_label, original_order)
plot1<- ggplot(for.plot, aes(x = c(1:dim(for.plot)[1]), y = power, color = ans_label)) +
  geom_point(size = 1, alpha = 0.7) +  # Adjust point sized
  ylim(0,1) +
  labs(title = "Overall power", x = "Individuals", y = "Power") +
  theme_minimal()
# ggsave("Overall_power_by_ancestry.png", plot = last_plot())


# plot homogeneous beta power along EUR percentage
for.plot.2 <- data.frame(power = power_homo_list, ans_label = df.true$ans_label)
ordered_indices <- order(df.true[,'EUR'])  # Ascending order of v
for.plot.2 <- for.plot.2[ordered_indices,]

plot2 <- ggplot(for.plot.2, aes(x = c(1:dim(for.plot.2)[1]), y = power, color = ans_label)) +
  geom_point(size = 1, alpha = 0.7) +  # Adjust point sized
  ylim(0,1) +
  labs(title = "Homogeneous effect power", x = "Individuals ordered by ascending EUR", y = "Power") +
  theme_minimal()
# ggsave("Homogeneous_power_by_ancestry.png", plot = last_plot(), height = 5, width = 8)

combined_plot <- plot1 + plot2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("combined_power_by_ancestry.png", plot = last_plot(), height = 5, width = 10)


library(ggplot2)
library(reshape2)

ancestry <- c(0.91, 0.14, 0.40, 0.69, 0.91, 0.20, 0.54, 0.33, 0.30, 0.29)
value1 <- c(1,2,3,4,5,6,7,8,9,10)
value2 <- c(1,2,3,4,5,6,7,8,9,10) * 0.1
df <- data.frame(ancestry = ancestry, value1 = value1, value2=value2)

ordered_indices <- order(df[,"ancestry"])
ordered_df <- df[ordered_indices, c("value1", "value2")]
rownames(ordered_df) <- NULL
long_matrix <- reshape2::melt(as.matrix(ordered_df))
colnames(long_matrix) <- c("Individual", "which.value", "Importance")

test_plot <- ggplot(long_matrix, aes(x = as.factor(Individual), y = which.value, fill = Importance)) +
  geom_tile() +
  geom_text(aes(label = round(Importance, 2)), color = "black", size = 2) +  # Add text labels
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal(base_size = 8) +
  guides(fill = guide_colorbar(keywidth = 0.3, keyheight = 2)) + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6),
    legend.position = "none",
    plot.title = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "MAF")


# Rachel's

ancestry <- c(0.91, 0.14, 0.40, 0.69, 0.91, 0.20, 0.54, 0.33, 0.30, 0.29)
value1 <- c(1,2,3,4,5,6,7,8,9,10)
value2 <- c(1,2,3,4,5,6,7,8,9,10) * 0.1
df <- data.frame(ancestry = ancestry, value1 = value1, value2=value2)

ordered_indices <- order(df[,"ancestry"])
ordered_df <- df[ordered_indices, c("value1", "value2")]

long_matrix <- reshape2::melt(as.matrix(ordered_df))
colnames(long_matrix) <- c("Individual", "which.value", "Importance")

long_matrix$Individual <- factor(long_matrix$Individual, levels = ordered_indices)

test_plot_2 <- ggplot(long_matrix, aes(x = Individual, y = which.value, fill = Importance)) +
  geom_tile() +
  geom_text(aes(label = round(Importance, 2)), color = "black", size = 2) +  # Add text labels
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal(base_size = 8) +
  guides(fill = guide_colorbar(keywidth = 0.3, keyheight = 2)) + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6),
    legend.position = "none",
    plot.title = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "MAF")

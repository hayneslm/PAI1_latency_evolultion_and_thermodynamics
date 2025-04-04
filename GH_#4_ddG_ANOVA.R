# Parses ddG values based on functional/non-functional in the emerical dataset. Makes box plots and performs ANOVA comparing the functional/nonfunctional variants.


create_boxplot_with_anova <- function(data, x_var, y_var, pos_var, title) {
  # Create groups based on the x_var condition
  data$group <- ifelse(data[[x_var]] > 0, ">0", "<=0")
  
  # Fit a linear model with ANOVA analysis
  formula <- as.formula(paste(y_var, "~ group +", pos_var))
  lm_result <- lm(formula, data = data)
  anova_result <- anova(lm_result)
  
  # Extract p-value for the group effect
  p_value <- anova_result$`Pr(>F)`[which(rownames(anova_result) == "group")]
  
  # Calculate means for each group
  group_means <- data %>%
    group_by(group) %>%
    summarize(mean_value = mean(get(y_var), na.rm = TRUE),
              count = n())
  
  mean_gt_0 <- group_means$mean_value[group_means$group == ">0"]
  mean_le_0 <- group_means$mean_value[group_means$group == "<=0"]
  count_gt_0 <- group_means$count[group_means$group == ">0"]
  count_le_0 <- group_means$count[group_means$group == "<=0"]
  
  # Handle potential NA p-values
  p_value_text <- if (!is.na(p_value)) {
    sprintf("p-value: %.2e", p_value)
  } else {
    "p-value: NA"
  }
  
  # Generate the boxplot with jittered points
  p <- ggplot(data, aes(x = factor(group), y = get(y_var), color = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ggtitle(title) +
    theme_classic() +
    labs(x = x_var, y = y_var, color = "Group") +
    annotate("text", x = 1, y = Inf, label = sprintf("Mean: %.2f\nN: %d", mean_le_0, count_le_0), hjust = 1.05, vjust = 1, size = 3.5) +
    annotate("text", x = 2, y = Inf, label = sprintf("Mean: %.2f\nN: %d", mean_gt_0, count_gt_0), hjust = 1.05, vjust = 1, size = 3.5) +
    annotate("text", x = 1.5, y = Inf, label = p_value_text, hjust = 0.5, vjust = 2, size = 3.5) +
    scale_color_manual(values = c("<=0" = "red", ">0" = "blue"))
  
  # Print ANOVA results to the console
  cat("\nANOVA results for", y_var, ":\n")
  print(anova_result)
  cat("\nCounts: ", y_var, " > 0: N = ", count_gt_0, ", ", y_var, " <= 0: N = ", count_le_0, "\n")
  
  return(p)
}

# Use the create_boxplot_with_anova function for analysis
plot_active_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "EvoEF_active", "Pos", "EvoEF Active 0h")
plot_latent_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "EvoEF_latent", "Pos", "EvoEF Latent 0h")
plot_AF_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "EvoEF_AF", "Pos", "EvoEF AF Conformation 0h")
plot_xfold_active_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "xfold_active", "Pos", "XFold Active Conformation 0h")
plot_xfold_latent_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "xfold_latent", "Pos", "XFold Latent Conformation 0h")
plot_xfold_AF_0h <- create_boxplot_with_anova(Compare_0h, "log2FoldChange", "xfold_AF", "Pos", "XFold AF Conformation 0h")

# Use the same function for 48h dataset
plot_active_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_active", "Pos", "EvoEF Active 48h")
plot_latent_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_latent", "Pos", "EvoEF Latent 48h")
plot_AF_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_AF", "Pos", "EvoEF AF Conformation 48h")
plot_xfold_active_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_active", "Pos", "FoldX Active Conformation 48h")
plot_xfold_latent_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_latent", "Pos", "FoldX Latent Conformation 48h")
plot_xfold_AF_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_AF", "Pos", "FoldX AF Conformation 48h")


# Arrange the 0h plots in a 3x2 grid
arranged_plots_0h <- ggarrange(
  plot_active_0h, plot_latent_0h, plot_AF_0h,
  plot_xfold_active_0h, plot_xfold_latent_0h, plot_xfold_AF_0h,
  ncol = 3, nrow = 2,
  common.legend = TRUE, legend = "bottom"
)

# Arrange the 48h plots in a 3x2 grid
arranged_plots_48h <- ggarrange(
  plot_active_48h, plot_latent_48h, plot_AF_48h,
  plot_xfold_active_48h, plot_xfold_latent_48h, plot_xfold_AF_48h,
  ncol = 3, nrow = 2,
  common.legend = TRUE, legend = "bottom"
)

# Print the arranged plots for 0h
print(arranged_plots_0h)

# Save the arranged plots for 0h to a file
ggsave("arranged_boxplots_0h.pdf", arranged_plots_0h, device = "pdf")

# Print the arranged plots for 48h
print(arranged_plots_48h)

# Save the arranged plots for 48h to a file
ggsave("arranged_boxplots_48h.pdf", arranged_plots_48h, device = "pdf")
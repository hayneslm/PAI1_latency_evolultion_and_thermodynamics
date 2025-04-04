# Parses VEP scores on the basis of functional/non-functional in the emperical dataset, performs ANOVA, and generates box plots.

library(ggplot2)
library(dplyr)
library(ggpubr)

create_boxplrescreate_boxplot_with_anova <- function(data, x_var, y_var, pos_var, title) {
  # Create groups based on the x_var condition
  data$group <- ifelse(data[[x_var]] > 0, ">0", "<=0")
  data$Pos <- as.factor(data$Pos)
  
  # Calculate means, standard errors, and counts for annotation
  calc_summary <- data %>%
    group_by(group) %>%
    summarise(
      mean = mean(get(y_var), na.rm = TRUE),
      n = n(),
      se = sd(get(y_var), na.rm = TRUE) / sqrt(n())
    )
  
  # Perform ANOVA
  anova_result <- aov(data[[y_var]] ~ data$group + data[[pos_var]])
  anova_summary <- summary(anova_result)
  p_value <- anova_summary[[1]]$`Pr(>F)`[1]
  
  # Generate the boxplot with jittered points
  p <- ggplot(data, aes(x = factor(group), y = get(y_var), color = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_text(data = calc_summary, aes(x = group, y = max(data[[y_var]], na.rm = TRUE), 
                                       label = sprintf("Mean: %.2f\nN: %d", mean, n)),
              hjust = 1.5, vjust = 1) +
    annotate("text", x = 1.5, y = max(data[[y_var]], na.rm = TRUE), label = sprintf("p-value: %.3f", p_value),
             hjust = 0.5, vjust = 4) +
    labs(title = title, x = x_var, y = y_var, color = "Group") +
    theme_bw() +
    scale_color_manual(values = c("<=0" = "red", ">0" = "blue"))
  
  # Create a summary data frame for ANOVA results
  anova_summary_df <- data.frame(
    Title = title,
    Group1_Mean = calc_summary$mean[calc_summary$group == "<=0"],
    Group2_Mean = calc_summary$mean[calc_summary$group == ">0"],
    p_value = p_value,
    Group1_N = calc_summary$n[calc_summary$group == "<=0"],
    Group2_N = calc_summary$n[calc_summary$group == ">0"]
  )
  
  return(list(plot = p, anova_summary = anova_summary_df))
}

# Generate boxplots and collect ANOVA summaries without printing or saving:
result_eve_0h <- create_boxplot_with_anova(EVE_score_0h, "log2FoldChange", "EVE_scores", "Pos", "EVE 0h")
result_am_0h <- create_boxplot_with_anova(AM_score_0h, "log2FoldChange", "AMscore", "Pos", "AM 0h")
result_cpt_0h <- create_boxplot_with_anova(CPT_score_0h, "log2FoldChange", "CPTscore", "Pos", "CPT 0h")
result_eve_48h <- create_boxplot_with_anova(EVE_score_48h, "log2FoldChange", "EVE_scores", "Pos", "EVE 48h")
result_am_48h <- create_boxplot_with_anova(AM_score_48h, "log2FoldChange", "AMscore", "Pos", "AM 48h")
result_cpt_48h <- create_boxplot_with_anova(CPT_score_48h, "log2FoldChange", "CPTscore", "Pos", "CPTscore 48h")

# Combine plots
plots <- list(
  result_eve_0h, result_am_0h, result_cpt_0h,
  result_eve_48h, result_am_48h, result_cpt_48h
)

# Arrange and save the plots together
arranged_plots <- ggarrange(plotlist = plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
print(arranged_plots)

# Save the composite plot to a file
ggsave("combined_boxplots.pdf", arranged_plots, width = 12, height = 12)

# Combine ANOVA summaries
anova_summaries <- bind_rows(
  result_eve_0h$plot_env$anova_result %>% mutate(source = "result_eve_0h"),
  result_am_0h$plot_env$anova_result %>% mutate(source = "result_am_0h"),
  result_cpt_0h$plot_env$anova_result %>% mutate(source = "result_cpt_0h"),
  result_eve_48h$plot_env$anova_result %>% mutate(source = "result_eve_48h"),
  result_am_48h$plot_env$anova_result %>% mutate(source = "result_am_48h"),
  result_cpt_48h$plot_env$anova_result %>% mutate(source = "result_cpt_48h")
)

# Save ANOVA summaries to a CSV file
write_csv(anova_summaries, "anova_summaries.csv")


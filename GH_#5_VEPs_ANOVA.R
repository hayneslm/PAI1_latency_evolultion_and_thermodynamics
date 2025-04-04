# Parses VEP scores on the basis of functional/non-functional in the emperical dataset, performs ANOVA, and generates box plots.

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(readr)
library(broom)
library(gridExtra)

EVE_score = read_csv("PAI1_human_eve.csv") %>% filter(mutations != "wt") %>% 
  separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(Pos = as.double(Pos)) %>% 
  mutate(Pos = Pos-23) %>% 
  mutate(Variant = paste(AA,Pos,Mut, sep = "")) %>% filter(Pos>0)

AM_score = read_tsv("PAI1_alphamissense.txt", col_names = F) %>% 
  select(!X1) %>% rename(mutations = X2, AMscore = X3, type = X4) %>% 
  separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(Pos = as.double(Pos)) %>% 
  mutate(Pos = Pos-23) %>% 
  mutate(Variant = paste(AA,Pos,Mut, sep = "")) %>% filter(Pos>0)

CPT_score = read_csv("PAI1_CPT_score.csv") %>%
  rename(mutations = mutant, CPTscore = CPT1_score) %>%  
  separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(Pos = as.double(Pos)) %>% 
  mutate(Pos = Pos-23) %>% 
  mutate(Variant = paste(AA,Pos,Mut, sep = "")) %>% filter(Pos>0) %>%
  filter(AA != Mut)


EVE_score_0h = full_join(Input_0h %>% mutate(Pos=as.numeric(Pos)),
                         EVE_score) %>% na.omit()

AM_score_0h = full_join(Input_0h %>% mutate(Pos=as.numeric(Pos)),
                        AM_score, by = c("Variant", "AA", "Pos", "Mut")) %>% 
  na.omit()

CPT_score_0h = full_join(Input_0h %>% mutate(Pos=as.numeric(Pos)),
                         CPT_score, by = c("Variant", "AA", "Pos", "Mut")) %>% 
  na.omit()


model <- lm(EVE_scores ~ log2FoldChange, data = EVE_score_0h)
summary(model)

ggplot(data = EVE_score_0h, aes(x = log2FoldChange, y = EVE_scores))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("EVE vs 0h data")+
  theme_classic()

model <- lm(AMscore ~ log2FoldChange, data = AM_score_0h)
summary(model)

ggplot(data = AM_score_0h, aes(x = log2FoldChange, y = AMscore))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("AM vs 0h data")+
  theme_classic()

model <- lm(CPTscore ~ log2FoldChange, data = CPT_score_0h)
summary(model)

ggplot(data = CPT_score_0h, aes(x = log2FoldChange, y = CPTscore))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("CPT vs 0h data")+
  theme_classic()
###########

EVE_score_48h = full_join(clean_48h %>% mutate(Pos=as.numeric(Pos)),
                          EVE_score) %>% na.omit()

AM_score_48h = full_join(clean_48h %>% mutate(Pos=as.numeric(Pos)),
                         AM_score, by = c("Variant", "AA", "Pos", "Mut")) %>% 
  na.omit()

CPT_score_48h = full_join(clean_48h %>% mutate(Pos=as.numeric(Pos)),
                          CPT_score, by = c("Variant", "AA", "Pos", "Mut")) %>% 
  na.omit()


model <- lm(EVE_scores ~ log2FoldChange, data = EVE_score_48h)
summary(model)

ggplot(data = EVE_score_48h, aes(x = log2FoldChange, y = EVE_scores))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("EVE vs 48h data")+
  theme_classic()

model <- lm(AMscore ~ log2FoldChange, data = AM_score_48h)
summary(model)

ggplot(data = AM_score_48h, aes(x = log2FoldChange, y = AMscore))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("AM vs 48h data")+
  theme_classic()

model <- lm(CPTscore ~ log2FoldChange, data = CPT_score_48h)
summary(model)

ggplot(data = CPT_score_48h, aes(x = log2FoldChange, y = CPTscore))+
  geom_jitter(color = 'grey60')+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  ylim(0,1)+
  ggtitle("CPT vs 48h data")+
  theme_classic()






# Function to fit the linear model and get a summary with R^2
summarize_model <- function(data, response, predictor, dataset_name) {
  model <- lm(reformulate(predictor, response = response), data = data)
  model_summary <- summary(model)
  data.frame(
    dataset = dataset_name,
    response_variable = response,
    r_squared = model_summary$r.squared,
    intercept_estimate = model_summary$coefficients[1, "Estimate"],
    predictor_estimate = model_summary$coefficients[2, "Estimate"],
    intercept_p_value = model_summary$coefficients[1, "Pr(>|t|)"],
    predictor_p_value = model_summary$coefficients[2, "Pr(>|t|)"]
  )
}

# # Running the linear models and capturing their summaries
model_summary_eve_0h <- summarize_model(EVE_score_0h, "EVE_scores", "log2FoldChange", "EVE_0h")
model_summary_am_0h <- summarize_model(AM_score_0h, "AMscore", "log2FoldChange", "AM_0h")
model_summary_cpt_0h <- summarize_model(CPT_score_0h, "CPTscore", "log2FoldChange", "CPT_0h")
model_summary_eve_48h <- summarize_model(EVE_score_48h, "EVE_scores", "log2FoldChange", "EVE_48h")
model_summary_am_48h <- summarize_model(AM_score_48h, "AMscore", "log2FoldChange", "AM_48h")
model_summary_cpt_48h <- summarize_model(CPT_score_48h, "CPTscore", "log2FoldChange", "CPT_48h")

# # Combining all summaries into one data frame
all_model_summaries <- bind_rows(model_summary_eve_0h, model_summary_am_0h,
                                 model_summary_cpt_0h,
                                 model_summary_eve_48h, model_summary_am_48h,
                                 model_summary_cpt_48h)
# 
# # Writing the summary to a CSV file
write_csv(all_model_summaries, "model_summaries.csv")
# 
# 


###########Make a Single Figure####################
p1 <- ggplot(data = EVE_score_0h, aes(x = log2FoldChange, y = EVE_scores)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("0h - EVE") +
  theme_classic()

p2 <- ggplot(data = AM_score_0h, aes(x = log2FoldChange, y = AMscore)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("0h - AM") +
  theme_classic()

p3 <- ggplot(data = CPT_score_0h, aes(x = log2FoldChange, y = CPTscore)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("0h - CPT") +
  theme_classic()

p4 <- ggplot(data = EVE_score_48h, aes(x = log2FoldChange, y = EVE_scores)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("48h - EVE") +
  theme_classic()

p5 <- ggplot(data = AM_score_48h, aes(x = log2FoldChange, y = AMscore)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("48h - AM") +
  theme_classic()

p6 <- ggplot(data = CPT_score_48h, aes(x = log2FoldChange, y = CPTscore)) +
  geom_jitter(color = 'grey60') +
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed") +
  ylim(0,1) +
  ggtitle("48h - CPT") +
  theme_classic()

# Open a PDF device to save the plots
pdf("VESPs_combined_plots.pdf", width = 15, height = 10)

# Combine plots
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

# Close the PDF device
dev.off()

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


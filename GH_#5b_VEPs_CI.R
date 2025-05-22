EVE_0h_split = EVE_score_0h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

EVE_0h_calc_summary <- EVE_0h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(EVE_scores),
    n = n(),
    se = sd(EVE_scores, na.rm = TRUE) / sqrt(n()))


EVE_0h_anova_result <- aov(EVE_0h_split$EVE_scores ~ EVE_0h_split$group + EVE_0h_split$Pos)
EVE_0h_anova_summary <- summary(EVE_0h_anova_result)
EVE_0h_anova_summary_p_value <- EVE_0h_anova_summary[[1]]$`Pr(>F)`[1]
confint(EVE_0h_anova_result)

t.test(EVE_0h_split$EVE_scores[EVE_0h_split$group == ">0"])
t.test(EVE_0h_split$EVE_scores[EVE_0h_split$group == "<=0"])


##
AM_48h_split = AM_score_48h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

AM_48h_calc_summary <- AM_48h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(AMscore),
    n = n(),
    se = sd(AMscore, na.rm = TRUE) / sqrt(n()))


AM_48h_anova_result <- aov(AM_48h_split$AMscore ~ AM_48h_split$group + AM_48h_split$Pos)
AM_48h_anova_summary <- summary(AM_48h_anova_result)
AM_48h_anova_summary_p_value <- AM_48h_anova_summary[[1]]$`Pr(>F)`[1]

t.test(AM_48h_split$AMscore[AM_48h_split$group == ">0"])
t.test(AM_48h_split$AMscore[AM_48h_split$group == "<=0"])

  
##
CPT_48h_split = CPT_score_48h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

CPT_48h_calc_summary <- CPT_48h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(CPTscore),
    n = n(),
    se = sd(CPTscore, na.rm = TRUE) / sqrt(n()))


CPT_48h_anova_result <- aov(CPT_48h_split$CPTscore ~ CPT_48h_split$group + CPT_48h_split$Pos)
CPT_48h_anova_summary <- summary(CPT_48h_anova_result) %>% print()
CPT_48h_anova_summary_p_value <- CPT_48h_anova_summary[[1]]$`Pr(>F)`[1]

t.test(CPT_48h_split$CPTscore[CPT_48h_split$group == ">0"])
t.test(CPT_48h_split$CPTscore[CPT_48h_split$group == "<=0"])

###########
EVE_48h_split = EVE_score_48h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

EVE_48h_calc_summary <- EVE_48h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(EVE_scores),
    n = n(),
    se = sd(EVE_scores, na.rm = TRUE) / sqrt(n()))


EVE_48h_anova_result <- aov(EVE_48h_split$EVE_scores ~ EVE_48h_split$group + EVE_48h_split$Pos)
EVE_48h_anova_summary <- summary(EVE_48h_anova_result) %>% print()
EVE_48h_anova_summary_p_value <- EVE_48h_anova_summary[[1]]$`Pr(>F)`[1]
confint(EVE_48h_anova_result)

t.test(EVE_48h_split$EVE_scores[EVE_48h_split$group == ">0"])
t.test(EVE_48h_split$EVE_scores[EVE_48h_split$group == "<=0"])


##
AM_48h_split = AM_score_48h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

AM_48h_calc_summary <- AM_48h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(AMscore),
    n = n(),
    se = sd(AMscore, na.rm = TRUE) / sqrt(n()))


AM_48h_anova_result <- aov(AM_48h_split$AMscore ~ AM_48h_split$group + AM_48h_split$Pos)
AM_48h_anova_summary <- summary(AM_48h_anova_result) %>% print()
AM_48h_anova_summary_p_value <- AM_48h_anova_summary[[1]]$`Pr(>F)`[1]

t.test(AM_48h_split$AMscore[AM_48h_split$group == ">0"])
t.test(AM_48h_split$AMscore[AM_48h_split$group == "<=0"])


##
CPT_48h_split = CPT_score_48h %>% mutate(
  group = ifelse(log2FoldChange > 0, ">0", "<=0"),
  Pos = as.factor(Pos))

CPT_48h_calc_summary <- CPT_48h_split %>%
  group_by(group) %>%
  summarize(
    mean = mean(CPTscore),
    n = n(),
    se = sd(CPTscore, na.rm = TRUE) / sqrt(n()))


CPT_48h_anova_result <- aov(CPT_48h_split$CPTscore ~ CPT_48h_split$group + CPT_48h_split$Pos)
CPT_48h_anova_summary <- summary(CPT_48h_anova_result) %>% print()
CPT_48h_anova_summary_p_value <- CPT_48h_anova_summary[[1]]$`Pr(>F)`[1]

t.test(CPT_48h_split$CPTscore[CPT_48h_split$group == ">0"])
t.test(CPT_48h_split$CPTscore[CPT_48h_split$group == "<=0"])


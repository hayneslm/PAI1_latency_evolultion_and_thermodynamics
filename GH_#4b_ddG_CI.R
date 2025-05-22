# EvoEF 0h active
  EvoEF_0hA_split = Compare_0h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  EvoEF_0hA_calc_summary <- EvoEF_0hA_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(EvoEF_active),
      n = n(),
      se = sd(EvoEF_active, na.rm = TRUE) / sqrt(n()))
  
  t.test(EvoEF_0hA_split$EvoEF_active[EvoEF_0hA_split$group == ">0"])
  t.test(EvoEF_0hA_split$EvoEF_active[EvoEF_0hA_split$group == "<=0"])

  # xfold 0h active
  xfold_0hA_split = Compare_0h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  xfold_0hA_calc_summary <- xfold_0hA_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(xfold_active),
      n = n(),
      se = sd(xfold_active, na.rm = TRUE) / sqrt(n()))
  
  t.test(xfold_0hA_split$xfold_active[xfold_0hA_split$group == ">0"])
  t.test(xfold_0hA_split$xfold_active[xfold_0hA_split$group == "<=0"])  
  
  
  
#####
  # EvoEF 0h latent
  EvoEF_0hL_split = Compare_0h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  EvoEF_0hL_calc_summary <- EvoEF_0hL_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(EvoEF_latent),
      n = n(),
      se = sd(EvoEF_latent, na.rm = TRUE) / sqrt(n()))
  
  t.test(EvoEF_0hL_split$EvoEF_latent[EvoEF_0hL_split$group == ">0"])
  t.test(EvoEF_0hL_split$EvoEF_latent[EvoEF_0hL_split$group == "<=0"])
  
  # xfold 0h latent
  xfold_0hL_split = Compare_0h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  xfold_0hL_calc_summary <- xfold_0hL_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(xfold_latent),
      n = n(),
      se = sd(xfold_latent, na.rm = TRUE) / sqrt(n()))
  
  t.test(xfold_0hL_split$xfold_latent[xfold_0hL_split$group == ">0"])
  t.test(xfold_0hL_split$xfold_latent[xfold_0hL_split$group == "<=0"])    
  
48h######################48h
  # EvoEF 48h active
  EvoEF_48hA_split = Compare_48h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  EvoEF_48hA_calc_summary <- EvoEF_48hA_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(EvoEF_active),
      n = n(),
      se = sd(EvoEF_active, na.rm = TRUE) / sqrt(n()))
  
  t.test(EvoEF_48hA_split$EvoEF_active[EvoEF_48hA_split$group == ">0"])
  t.test(EvoEF_48hA_split$EvoEF_active[EvoEF_48hA_split$group == "<=0"])
  
  # xfold 48h active
  xfold_48hA_split = Compare_48h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  xfold_48hA_calc_summary <- xfold_48hA_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(xfold_active),
      n = n(),
      se = sd(xfold_active, na.rm = TRUE) / sqrt(n()))
  
  t.test(xfold_48hA_split$xfold_active[xfold_48hA_split$group == ">0"])
  t.test(xfold_48hA_split$xfold_active[xfold_48hA_split$group == "<=0"])  
  
  
  
  #####
  # EvoEF 48h latent
  EvoEF_48hL_split = Compare_48h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  EvoEF_48hL_calc_summary <- EvoEF_48hL_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(EvoEF_latent),
      n = n(),
      se = sd(EvoEF_latent, na.rm = TRUE) / sqrt(n()))
  
  t.test(EvoEF_48hL_split$EvoEF_latent[EvoEF_48hL_split$group == ">0"])
  t.test(EvoEF_48hL_split$EvoEF_latent[EvoEF_48hL_split$group == "<=0"])
  
  # xfold 48h latent
  xfold_48hL_split = Compare_48h %>% mutate(
    group = ifelse(log2FoldChange > 0, ">0", "<=0"),
    Pos = as.factor(Pos))
  
  xfold_48hL_calc_summary <- xfold_48hL_split %>%
    group_by(group) %>%
    summarize(
      mean = mean(xfold_latent),
      n = n(),
      se = sd(xfold_latent, na.rm = TRUE) / sqrt(n()))
  
  t.test(xfold_48hL_split$xfold_latent[xfold_48hL_split$group == ">0"])
  t.test(xfold_48hL_split$xfold_latent[xfold_48hL_split$group == "<=0"])    
  
  
  
  
    
  # Use the same function for 48h dataset
  plot_active_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_active", "Pos", "EvoEF Active 48h")
  plot_latent_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_latent", "Pos", "EvoEF Latent 48h")
  plot_AF_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "EvoEF_AF", "Pos", "EvoEF AF Conformation 48h")
  plot_xfold_active_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_active", "Pos", "FoldX Active Conformation 48h")
  plot_xfold_latent_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_latent", "Pos", "FoldX Latent Conformation 48h")
  plot_xfold_AF_48h <- create_boxplot_with_anova(Compare_48h, "log2FoldChange", "xfold_AF", "Pos", "FoldX AF Conformation 48h")
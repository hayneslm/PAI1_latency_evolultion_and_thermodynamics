# Compare conservation scores to emperically detemined DMS data

library(ggplot2)
library(dplyr)
library(readr)

setwd("/Users/hayneslm/Documents/R_projects/Binning_2024-06/")

###Compare Conservation Scores###
conservScore = NULL
conservScore = read_tsv("data_scores.txt") %>% rename(Pos = Position)
all_pos = c(1:379) %>% as.numeric() %>% as_tibble() %>% rename(Pos = value) %>% 
  mutate(score = "score.x")

clean_48h = clean_48h #%>% filter(padj < 0.1)

h48_count_total = clean_48h %>% group_by(Pos) %>% 
  summarize(count = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(count = as.numeric(count))

h48_counts_enriched = clean_48h %>% filter(log2FoldChange>0) %>% 
  group_by(Pos) %>% 
  summarize(countEC = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(countEC = as.numeric(countEC))
h48_norm_count = NULL
h48_norm_count = full_join(h48_count_total, 
                           h48_counts_enriched, by = "Pos") %>% 
  mutate(h48_norm_count = countEC/count)

h48_counts_depleted = clean_48h %>% filter(log2FoldChange<=0) %>%
  group_by(Pos) %>%
  summarize(countEC = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(countEC = as.numeric(countEC))
h48_norm_count_dep = NULL
h48_norm_count_dep = full_join(h48_count_total,
                           h48_counts_depleted, by = "Pos") %>%
  mutate(h48_norm_count_dep = countEC/count)


h0_count_total = Input_0h_filtered %>% group_by(Pos) %>% 
  summarize(count = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(count = as.numeric(count))

h0_counts_enriched = Input_0h_filtered %>% filter(log2FoldChange>0) %>% 
  group_by(Pos) %>% 
  summarize(countEC = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(countEC = as.numeric(countEC))
h0_norm_count = NULL
h0_norm_count = full_join(h0_count_total, 
                           h0_counts_enriched, by = "Pos") %>% 
  mutate(h0_norm_count = countEC/count)

h0_counts_depleted = Input_0h_filtered %>% filter(log2FoldChange<=0) %>% 
  group_by(Pos) %>% 
  summarize(countEC = n()) %>% mutate(Pos = as.numeric(Pos)) %>%
  mutate(countEC = as.numeric(countEC))
h0_norm_count_dep = NULL
h0_norm_count_dep = full_join(h0_count_total, 
                          h0_counts_depleted, by = "Pos") %>% 
  mutate(h0_norm_count_dep = countEC/count)


conservScore = full_join(conservScore,h48_norm_count, by = "Pos")
conservScore = full_join(conservScore, h0_norm_count, by = "Pos")
conservScore = full_join(conservScore,h48_norm_count_dep, by = "Pos")
conservScore = full_join(conservScore, h0_norm_count_dep, by = "Pos")

# h48_conservScore = conservScore %>% select(h48_norm_count, score) %>%
#   na.omit() %>%
#   mutate(h48_quartile = cut_number(h48_norm_count, n = 4))

# ggplot(h48_conservScore, aes(x = h48_quartile, y = score,
#                              color = h48_quartile))+
#   geom_violin()+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   theme_classic()+
#   theme(legend.position = "none")+geom_boxplot(width=0.1)+
#   ggtitle("48h enriched")

# write_tsv(h48_conservScore,"h48_conservScore.xls")

# h0_conservScore = conservScore %>% select(h0_norm_count, score) %>%
#   na.omit() %>%
#   mutate(h0_quartile = cut_number(h0_norm_count, n = 4))
# 
# ggplot(h0_conservScore, aes(x = h0_quartile, y = score,
#                              color = h0_quartile))+
#   geom_violin()+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   theme_classic()+
#   theme(legend.position = "none")+geom_boxplot(width=0.1)+
#   ggtitle("0h enriched")
# 
# write_tsv(h0_conservScore,"h0_conservScore.xls")

# h0_conservScore_dep = conservScore %>% select(h0_norm_count_dep, score) %>%
  # na.omit() %>%
  # mutate(h0_quartile = cut_number(h0_norm_count_dep, n = 4))

# ggplot(h0_conservScore_dep, aes(x = h0_quartile, y = score,
#                             color = h0_quartile))+
#   geom_violin()+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   theme_classic()+
#   theme(legend.position = "none")+geom_boxplot(width=0.1)+
#   ggtitle("0h depleted")
# 
# h48_conservScore_dep = conservScore %>% select(h48_norm_count_dep, score) %>%
#   na.omit() %>%
#   mutate(h48_quartile = cut_number(h48_norm_count_dep, n = 4))
# 
# ggplot(h48_conservScore_dep, aes(x = h48_quartile, y = score,
#                             color = h48_quartile))+
#   geom_violin()+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   theme_classic()+
#   theme(legend.position = "none")+geom_boxplot(width=0.1)+
#   ggtitle("48h depleted")


# Compare h48 to h0 colored with dN
Compare_48h_0h_conServ = conservScore %>% 
  select(Pos, h0_norm_count, h48_norm_count, score) %>% 
  na.omit() %>% 
  mutate(quad = if_else(h0_norm_count>median(Compare_48h_0h_conServ$h0_norm_count) & h48_norm_count>median(Compare_48h_0h_conServ$h48_norm_count),"quad_1",
                        if_else(h0_norm_count<=median(Compare_48h_0h_conServ$h0_norm_count) & h48_norm_count>median(Compare_48h_0h_conServ$h48_norm_count), "quad_2",
                                if_else(h0_norm_count<=median(Compare_48h_0h_conServ$h0_norm_count) & h48_norm_count<=median(Compare_48h_0h_conServ$h48_norm_count), "quad_3", "quad_4"))))

write_tsv(Compare_48h_0h_conServ, "Compare_48h_0h_conServ.xls")

# Compare_48h_0h_conServ_plusVtn = read_tsv("Compare_48h_0h_conServ_plusVtn.txt") %>% rename(VtnCat = "...6") #%>% mutate(VtnCat = as.factor(VtnCat))
                       
# ggplot(Compare_48h_0h_conServ_plusVtn, 
#        aes(x = h0_norm_count, y = h48_norm_count, size = 1, shape = VtnCat))+
#   geom_jitter()+scale_color_gradientn(colours = topo.colors(6))+
#   theme_classic()+
#   geom_hline(yintercept = 0.5, linetype = "dotted")+
#   geom_vline(xintercept = 0.5, linetype = "dotted")


# ggplot(Compare_48h_0h_conServ_plusVtn, 
#        aes(x = h0_norm_count, y = h48_norm_count, size = 1, shape = VtnCat))+
#   geom_jitter()+
#   scale_color_gradientn(colours = topo.colors(6))+
#   theme_classic()+
#   geom_hline(yintercept = 0.5, linetype = "dotted")+
#   geom_vline(xintercept = 0.5, linetype = "dotted")




# ggplot(Compare_48h_0h_conServ_plusVtn, # %>% filter(VtnCat == "C"),
#        aes(x = h0_norm_count, y = h48_norm_count, 
#            color = score, size = VtnCat, shape = VtnCat))+
#   geom_jitter()+scale_color_gradientn(colours = topo.colors(6))+
#   theme_classic()+
#   geom_hline(yintercept = 0.5, linetype = "dotted")+
#   geom_vline(xintercept = 0.5, linetype = "dotted")

# 
ggplot(Compare_48h_0h_conServ, aes(x = quad, y = score, color = quad))+
  geom_violin()+
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, alpha = 0.5)+
  theme_classic()+
  geom_boxplot(width = 0.1)

ggsave("quad_compar_violin.pdf", device = "pdf")

# 
compare = aov(lm(score~quad, data = Compare_48h_0h_conServ))
summary(compare)
TukeyHSD(compare)

summary_stats <- Compare_48h_0h_conServ %>%
  group_by(quad) %>%
  summarise(
    count = n(),
    mean_score = mean(score, na.rm = TRUE),
    range_score = paste(range(score, na.rm = TRUE), collapse = " to "), # Format range as a string
    median_score = median(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE)
  )

print(summary_stats)
  
  
# ggplot(Compare_48h_0h_conServ_plusVtn, aes(x = VtnCat, y = score, color = VtnCat))+
#   geom_violin()+
#   geom_jitter(shape = 16, position = position_jitter(0.2))+
#   theme_classic()+
#   geom_boxplot(width = 0.1)
# 
# 
# Compare_48h_0h_conServ_hl = Compare_48h_0h_conServ %>% 
#   mutate(h0_hl = if_else(quad %in% c("quad_1","quad_4"),"high","low")) %>% 
#   mutate(h48_hl = if_else(quad %in% c("quad_1","quad_2"),"high","low"))
# 
# compare = aov(lm(score~h0_hl+h48_hl, data = Compare_48h_0h_conServ_hl))
# summary(compare)
# TukeyHSD(compare)
#  

compare = aov(lm(score~quad, data = Compare_48h_0h_conServ))
summary(compare)
TukeyHSD(compare)

# model <- lm(h48_norm_count ~ h0_norm_count, data = Compare_48h_0h_conServ)
# summary(model)


# Reorder the dataset by score in ascending order (smallest to largest)
Compare_48h_0h_conServ_ordered <- Compare_48h_0h_conServ %>%
  arrange(score)

# Create the plot with the reordered data
ggplot(Compare_48h_0h_conServ_ordered, 
       aes(x = h0_norm_count, y = h48_norm_count, color = score)) +
  geom_jitter(size = 6) +
  scale_color_gradientn(colours = topo.colors(10)) +
  theme_classic() +
  geom_hline(yintercept = median(Compare_48h_0h_conServ$h48_norm_count), linetype = "dotted") +
  geom_vline(xintercept = median(Compare_48h_0h_conServ$h0_norm_count), linetype = "dotted")
ggsave("Compare_quads.pdf", device = "pdf")

write_tsv(Compare_48h_0h_conServ, "Compare_48h_0h_conServ.xls")

# This is just looking at the regression in quads 1 and 2 but I don't think it is tellin us much of anything
# 
# test = Compare_48h_0h_conServ %>% filter(quad %in% c("quad_2","quad_3"))
# 
# plot(test$h48_norm_count, test$score)
# model <- lm(score ~ h48_norm_count, data = test)
# summary(model)
# 
# ggplot(test, aes(x = h48_norm_count, y = score, color = Pos)) +
#   geom_jitter(size = 4) +
#   geom_smooth(method = "lm", se = FALSE, color = 'black', linetype = "dashed") +
#   scale_color_gradientn(colours = terrain.colors(10)) +
#   theme_classic()

# 

compare = aov(lm(score~h0_norm_count*h48_norm_count, data = Compare_48h_0h_conServ))
summary(compare)
TukeyHSD(compare)



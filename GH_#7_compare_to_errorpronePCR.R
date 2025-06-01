library(readr)

SR_0h_data = read_tsv("WTbg_0h_screen.txt") %>% rename(Variant = mutation)
JBC_48h_data = read_csv("jbc2022_48h_data_filtered.csv") %>% rename(Variant = mutation)

# Compare the 0h data between the two data sets
Compare2old_0h = full_join(Input_0h_filtered, SR_0h_data, by = "Variant") %>% na.omit()

model <- lm(log2FoldChange.y ~ log2FoldChange.x, data = Compare2old_0h)
summary(model)

ggplot(Compare2old_0h, aes(x = log2FoldChange.x, y = log2FoldChange.y))+
  geom_jitter(color = "grey60")+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  theme_classic()

ggsave("Compare2old_0h.pdf", device = "pdf")

# Compare the 48h data between the two data sets
Compare2old_48h = full_join(Input_48h_filtered, JBC_48h_data, by = "Variant") %>% na.omit()

model <- lm(log2FoldChange.y ~ log2FoldChange.x, data = Compare2old_48h)
summary(model)

ggplot(Compare2old_48h, aes(x = log2FoldChange.x, y = log2FoldChange.y))+
  geom_jitter(color = "grey60")+
  geom_smooth(method = "lm", se = F, color = 'black', linetype = "dashed")+
  theme_classic()

ggsave("Compare2old_48h.pdf", device = "pdf")

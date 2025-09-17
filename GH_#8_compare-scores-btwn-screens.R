# This is a posthoc analysis to compare scores at 0h and 48h screens. "Clean" refers to those variants that were identified functional in the 0h screen and "Dirty" refers to those variants that were non-functional at 0h and therefore not carried through to the 48h screen analysis.

dirty_48h = anti_join(Input_48h_filtered, clean_48h, by = "Variant") %>% 
  rbind(.,Input_48h %>% filter(type == "wt"))

clean_48h_b = clean_48h %>% rbind(.,Input_48h %>% filter(type == "wt"))

library(ggplot2)

ggplot(filtered, aes(x = log2FoldChange, fill = type2)) +
  geom_histogram(binwidth = .1, position = "identity", alpha = 0.5) +
  scale_fill_manual(
    values = c("Dirty" = "dodgerblue", "Clean" = "forestgreen"),
    name = "Sample Type"
  ) +
  labs(
    title = "Distribution of log2FoldChange by Sample Type (No WT, No Score)",
    x = "log2FoldChange",
    y = "Count"
  ) +
  theme_classic()

ggsave("distribution_clean_dirty.pdf", device = "pdf")

shapiro.test(filtered$log2FoldChange[filtered$type2 == "Dirty"])
shapiro.test(filtered$log2FoldChange[filtered$type2 == "Clean"])

wilcox.test(log2FoldChange ~ type2, data = filtered)

aggregate(log2FoldChange ~ type2, data = filtered, function(x) c(n=length(x), mean=mean(x), median=median(x)))





dirty_48h_list = dirty_48h %>% select(Variant)

t0_values = semi_join(Input_0h, dirty_48h_list, by = "Variant")

t48_values = semi_join(Input_48h, dirty_48h_list, by = "Variant")

values_combined = full_join(t0_values, t48_values, by = "Variant") %>% 
  rename(log2FoldChange_0h = log2FoldChange.x, log2FoldChange_48h = log2FoldChange.y, type = type.x) #%>% 
  filter(Variant != "I91I")
  
# Plot 0h vs 48h, color by whether 'type' is wt or not
ggplot(values_combined, aes(x = log2FoldChange_0h, y = log2FoldChange_48h)) +
  geom_point(aes(color = type == "wt"), size = 2, alpha = 0.8) +
  scale_color_manual(
    name = "is WT?",
    values = c("TRUE" = "gold", "FALSE" = "gray60"),
    labels = c("FALSE" = "Other", "TRUE" = "WT")
  ) +
  labs(
    title = "log2FoldChange at 0h vs 48h",
    x = "log2FoldChange 0h",
    y = "log2FoldChange 48h"
  ) +
  theme_minimal()

library(ggrepel)

ggplot(values_combined, aes(x = log2FoldChange_0h, y = log2FoldChange_48h)) +
  geom_point(aes(color = type == "wt"), size = 2, alpha = 0.8) +
  geom_text_repel(
    data = subset(values_combined, type == "wt"),
    aes(label = Variant),  # assuming 'Variant' is a column in your frame
    color = "goldenrod", size = 3
  ) +
  scale_color_manual(
    name = "is WT?",
    values = c("TRUE" = "gold", "FALSE" = "gray60"),
    labels = c("FALSE" = "Other", "TRUE" = "WT")
  ) +
  labs(
    title = "log2FoldChange at 0h vs 48h",
    x = "log2FoldChange 0h",
    y = "log2FoldChange 48h"
  ) +
  theme_minimal()





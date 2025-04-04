# This code performs an initial analysis of the data generated in GH_#1_DeSeq.R, including MA plots and variant heatmaps.


# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("pheatmap")
# install.packages("reshape2")
library(pheatmap)
library(dplyr)
library(tidyverse)
library(readr)
library(scales)
library(reshape2)

setwd("~/Documents/R_projects/Binning_2024-06")

# Reimport DeSeq2 results
Input_0h = read_tsv("condition_h0_vs_Input_total.txt") %>% 
  na.omit() %>%
  mutate(score = if_else(padj<0.1, "pass", "fail"))%>% 
  mutate(temp = Variant) %>%   
  separate(temp, into = c("AA","Rest"),sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into = c("Pos", "Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(type = if_else(AA==Mut, "wt",
                        if_else(Mut=="X", "nonsense",
                                if_else(Mut=="B", "amber","missense"))))

Input_48h = read_tsv("condition_h48_vs_Input_total.txt") %>% 
  na.omit() %>%
  mutate(score = if_else(padj<0.1, "pass", "fail"))%>% 
  mutate(temp = Variant) %>%   
  separate(temp, into = c("AA","Rest"),sep = "(?<=[A-Z])(?=[0-9])") %>% 
  separate(Rest, into = c("Pos", "Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(type = if_else(AA==Mut, "wt",
                        if_else(Mut=="X", "nonsense",
                                if_else(Mut=="B", "amber","missense"))))


# MA Plots

x_breaks <- c(1e1, 1e3, 1e5, 1e7)

ggplot()+
  scale_x_log10(breaks = x_breaks)+
  geom_hline(yintercept = 0, color = "black") +
  xlab("Base Mean Score")+
  ylab(expression(bold(Log[2] ~ "-Fold Change")))+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_0h, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_0h, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_0h, type == "wt"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_0h, type == "amber"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_vline(xintercept = 50, color = "black", linetype = "dashed")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14), # Change the legend title font
    legend.text = element_text(size = 12))+
  guides(shape = guide_none()) + # This line will remove the legend for "score"
  labs(color = "Mutation Type") # Rename 'type' legend title to 'Mutation Type'
  
ggsave("0h_MA.pdf", device = pdf)


ggplot()+
  scale_x_log10(breaks = x_breaks)+
  geom_hline(yintercept = 0, color = "black") +
  xlab("Base Mean Score")+
  ylab(expression(bold(Log[2] ~ "-Fold Change")))+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_48h, type == "missense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_48h, type == "nonsense"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_48h, type == "wt"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_jitter(data = subset(Input_48h, type == "amber"), 
              aes(baseMean, log2FoldChange, color = type, shape = score), 
              stroke = 0.5, size = 1.5)+
  geom_vline(xintercept = 50, color = "black", linetype = "dashed")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14), # Change the legend title font
    legend.text = element_text(size = 12))+
  guides(shape = guide_none()) + # This line will remove the legend for "score"
  labs(color = "Mutation Type") # Rename 'type' legend title to 'Mutation Type'

ggsave("48h_MA.pdf", device = pdf)

#Set basemean cuttoff to 50
Input_0h_filtered = Input_0h %>% filter(#score == 'pass', 
                                        baseMean > 50,
                                          type == "missense")
Input_48h_filtered = Input_48h %>% filter(#score == 'pass', 
                                          baseMean > 50,
                                          type == "missense")

# Make Heatmaps
# Currently with lots of filters but they look better if you just let the data be the data. (took out filters)
#Input information for generating heatmaps
aa20 =  "AVLIMFYWRKHDESTNQGCP"
aa20 = unlist(strsplit(aa20, split = ""))

#use WT PAI-1 sequence
PAI1seq = "VHHPPSYVAHLASDFGVRVFQQVAQASKDRNVVFSPYGVASVLAMLQLTTGGETQQQIQAAMGFKIDDKGMAPALRHLYKELMGPWNKDEISTTDAIFVQRDLKLVQGFMPHFFRLFRSTVKQVDFSEVERARFIINDWVKTHTKGMISNLLGKGAVDQLTRLVLVNALYFNGQWKTPFPDSSTHRRLFHKSDGSTVSVPMMAQTNKFNYTEFTTPDGHYYDILELPYHGDTLSMFIAAPYEKEVPLSALTNILSAQLISHWKGNMTRLPRLLVLPKFSLETEVDLRKPLENLGMTDMFRQFQADFTSLSDQEPLHVAQALQKVKIEVNESGTVASSSTAVIVSARMAPEEIIMDRPFLFVVRHNPTGTVLFMGQVMEP"
Pos = c(1:nchar(PAI1seq))
AA = unlist(strsplit(PAI1seq, split=""))
WT_coded = cbind(Pos, AA, replicate(nchar(PAI1seq),10)) %>% as.data.frame() %>% 
  rename(value = V3) %>% mutate(value = as.numeric(value)) %>%
  rename(Mut = AA) %>% 
  mutate(padj = NA) %>% 
  as_tibble()

# Heatmap for 0h selection
# Input_0h_HM = Input_0h %>% filter(baseMean>=50)
Input_0h_HM = Input_0h_filtered
HM_data = Input_0h_HM %>% select(Pos,Mut,log2FoldChange,padj) %>% 
  rename(value = log2FoldChange)
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
# levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0, 
                                          max(HM_data$value), 10), 
                                        to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=30),
        axis.title = element_text(size = 36))
  # ggtitle("0h selection")


ggsave("0h_fingerprint.png", device = "png", 
       width = 44, height = 8.5, units = "in")

# Heatmap for 48h selection
# Input_48h_HM = Input_48h %>% filter(baseMean>=50)
Input_48h_HM = Input_48h_filtered
HM_data = Input_48h_HM %>% select(Pos,Mut,log2FoldChange,padj) %>% 
  rename(value = log2FoldChange)
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
# levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0, 
                                          max(HM_data$value), 10), 
                                        to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=30),
        axis.title = element_text(size = 36))
# ggtitle("0h selection")

ggsave("48h_fingerprint.png", device = "png", 
       width = 44, height = 8.5, units = "in")


# Make a 48h selected (clean_48h) that doesn't include non-functional at 0h
# nf_0h = Input_0h %>% filter(log2FoldChange <= 0) %>% select(Variant)
# clean_48h = anti_join(Input_48h, nf_0h, by = "Variant")

func_0h = Input_0h_filtered %>% filter(log2FoldChange > 0) %>% select(Variant)
clean_48h = right_join(Input_48h_filtered, func_0h) %>% na.omit() %>%
  filter(type == "missense")

wt_48h = Input_48h %>% filter(type == "wt")

# clean_48h_WTnorm = full_join(clean_48h, wt_48h, by = c("AA","Pos")) %>%
#   na.omit() %>% mutate(log2FoldChange = log2FoldChange.x-log2FoldChange.y)
# 
# #heatmap
# Input_48hClean_HM = clean_48h_WTnorm #%>% filter(baseMean>=50)
# HM_data = Input_48hClean_HM %>% select(Pos,Mut.x,log2FoldChange,padj.x) %>%
#   rename(value = log2FoldChange, Mut = Mut.x, padj = padj.x)
# HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>%
#   filter(Pos <= 379)
# HMlist     <- HMlist[order(HMlist$Pos),]
# HMlist$Pos <- as.factor(HMlist$Pos)
# # levels(HMlist$Pos)
# 
# ggplot(HMlist, aes(x = Pos, y = Mut)) +
#   geom_tile(aes(fill=value))+
#   scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
#   scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
#   scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
#                        values = rescale(c(min(HM_data$value), 0,
#                                           max(HM_data$value), 10),
#                                         to = c(0,1)),
#                        na.value = "ivory2") +
#   theme(panel.background = element_rect(fill = 'white'),
#         #legend.position = "none",
#         axis.text=element_text(size=24))+
#   ggtitle("48h selection clean")
# 
# ggsave("48h_clean_fingerprint.png", device = "png",
#        width = 44, height = 8.5, units = "in")

#heatmap
Input_48hClean_HM = clean_48h #%>% filter(baseMean>=50)
HM_data = Input_48hClean_HM %>% select(Pos,Mut,log2FoldChange,padj) %>%
  rename(value = log2FoldChange)
HMlist = rbind(HM_data, WT_coded) %>% mutate(Pos = as.numeric(Pos)) %>%
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)
# levels(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(colors = c("red", "snow2", "blue", "grey30"),
                       values = rescale(c(min(HM_data$value), 0,
                                          max(HM_data$value), 10),
                                        to = c(0,1)),
                       na.value = "ivory2") +
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "none",
        axis.text=element_text(size=30),
        axis.title = element_text(size = 36))
# ggtitle("0h selection")

ggsave("48h_clean_fingerprint.png", device = "png",
       width = 44, height = 8.5, units = "in")



# This is for a heatmap of the variants present in the library
count_data = all_amp_data_processed %>% 
  rename_with(~c("h48","h48","Input","h0","Input",
                 "h48","Input","h0","h0", "AA","Pos","Mut") %>% make.unique()) %>% 
  filter(AA != Mut, Mut != "B", Mut != "X")


count_data_Input = count_data %>% select(!contains(c("h48","h0"))) %>% 
  mutate(Input_sum = Input+Input.1+Input.2, 
         Input_mean = (Input+Input.1+Input.2)/3)

sum_sum = sum(count_data_Input$Input_sum)
sum_mean = sum(count_data_Input$Input_mean)

count_data_fraction = count_data_Input %>% mutate(sum_frac = Input_sum/sum_sum,
                                                  mean_frac = Input_mean/sum_mean) %>%
  mutate(mean_percent = mean_frac*100)
HM_data = count_data_fraction %>% select(Pos,Mut,mean_percent) %>% 
  rename(value = mean_percent)
HMlist = rbind(HM_data, WT_coded %>% select(!padj) %>% mutate(value = case_when(
  value == 10 ~ NA_real_, TRUE~value))) %>% mutate(Pos = as.numeric(Pos)) %>% 
  filter(Pos <= 379)
HMlist     <- HMlist[order(HMlist$Pos),]
HMlist$Pos <- as.factor(HMlist$Pos)

ggplot(HMlist, aes(x = Pos, y = Mut)) +
  geom_tile(aes(fill=value))+
  scale_x_discrete(name = "Amino Acid Position",  breaks=c(seq(0,379,50)))+
  scale_y_discrete(name = "Amino Acid", limit = rev(aa20))+
  scale_fill_gradientn(
    colors = c("white","grey100", "grey80","grey60","grey40", "grey20", "black"),
    values = rescale(c(min(HMlist$value, na.rm=TRUE), 
                       median(HMlist$value, na.rm=TRUE), 
                       max(HMlist$value, na.rm=TRUE))),
    na.value = "royalblue",
    name = "Percent of Library")+
  theme(panel.background = element_rect(fill = 'white'),
        legend.position = "right",
        axis.text=element_text(size=30),
        axis.title = element_text(size = 36))

ggsave("Library_distribution.pdf", device = "pdf", 
       width = 44, height = 8.5, units = "in")


# 

# clean_48h_sum = clean_48h %>% group_by(Pos) %>% 
#   summarize(posSum = sum(log2FoldChange) ) %>% 
#   mutate(Pos = as.numeric(Pos))
# 
# write_tsv(clean_48h_sum, "clean_48h_sum.xls")
# 
# ggplot(data = clean_48h_sum, aes(Pos, posSum)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "loess", se = FALSE, color = "blue", na.rm = TRUE)+
#   theme_classic()
# 
# ggplot(data = clean_48h_sum, aes(Pos, posSum)) +
#   geom_point(size = 3) +
#   geom_smooth(
#     method = "loess",
#     se = FALSE,
#     color = "blue",
#     na.rm = TRUE,
#     method.args = list(span = 0.1, degree = 2, surface = "direct", family = "symmetric")
#   ) +
#   theme_classic()


# # Calculate mean log2FoldChange per Pos
# clean_48h_mean <- clean_48h %>%
#   group_by(Pos) %>%
#   summarize(
#     posSum = sum(log2FoldChange),
#     count = n()
#   ) %>%
#   mutate(
#     posMean = posSum / count,
#     Pos = as.numeric(Pos)
#   )
# 
# # Write the results to a file
# write_tsv(clean_48h_mean, "clean_48h_mean.xls")
# 
# # Plot posMean vs Pos
# ggplot(data = clean_48h_mean, aes(Pos, posMean)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "loess", se = FALSE, color = "blue", na.rm = TRUE) +
#   theme_classic()
# 
# # Plot with additional smoothing arguments
# ggplot(data = clean_48h_mean, aes(Pos, posMean)) +
#   geom_point(size = 3) +
#   geom_smooth(
#     method = "loess",
#     se = FALSE,
#     color = "blue",
#     na.rm = TRUE,
#     method.args = list(span = 0.1, degree = 2, surface = "direct", family = "symmetric")
#   ) +
#   theme_classic()

# Calculate mean log2FoldChange per Pos, count of values, 
# count of values greater than 0, and the ratio of values greater than 0
clean_48h_summary <- clean_48h %>%
  group_by(Pos) %>%
  summarize(
    posSum = sum(log2FoldChange),
    count = n(),
    gt_zero_count = sum(log2FoldChange > 0)
  ) %>%
  mutate(
    posMean = posSum / count,
    Pos = as.numeric(Pos),
    ratio_gt_zero = gt_zero_count / count
  )

# Write the results to a file
write_tsv(clean_48h_summary, "clean_48h_summary.xls")

# Plot posMean vs Pos with a smooth line
# ggplot(data = clean_48h_summary, aes(Pos, posMean)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "loess", se = FALSE, color = "blue", na.rm = TRUE) 
#   theme_classic()

# # Plot posMean vs Pos with additional smoothing arguments
# ggplot(data = clean_48h_summary, aes(Pos, posMean)) +
#   geom_point(size = 3) +
#   geom_smooth(
#     method = "loess",
#     se = FALSE,
#     color = "blue",
#     na.rm = TRUE,
#     method.args = list(span = 0.1, degree = 2, surface = "direct", family = "symmetric")
#   ) +
#   theme_classic()

# # Plot the ratio of log2FoldChange > 0 vs Pos with a smooth line
# ggplot(data = clean_48h_summary, aes(Pos, ratio_gt_zero)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "loess", se = FALSE, color = "blue", na.rm = TRUE) +
#   theme_classic()
# 
# # Plot the ratio of log2FoldChange > 0 vs Pos with additional smoothing arguments
# ggplot(data = clean_48h_summary, aes(Pos, ratio_gt_zero)) +
#   geom_point(size = 3) +
#   geom_smooth(
#     method = "loess",
#     se = FALSE,
#     color = "blue",
#     na.rm = TRUE,
#     method.args = list(span = 0.1, degree = 2, surface = "direct", family = "symmetric")
#   ) +
#   theme_classic()


# This code install necessary packages and libraries and performs a DESeq2 analysis of the variant data as counts of each amino acid substitution in the PAI-1 input and 0h/48h selected variant libraries.


# install.packages("tidyverse");
# install.packages('dplyr')
# install.packages("reshape2");
# install.packages('factoextra')
# install.packages("ggpubr")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# install.packages('ggExtra')
# install.packages('RColorBrewer')
# install.packages("plot3D")
# install.packages("pheatmap)


library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(readr)
library(DESeq2)
library(ggpubr)
library(stats)
library(MASS)
library(plot3D)
library(tidyverse)
library(dplyr)

setwd("~/Documents/R_projects/Binning_2024-06/")
ampcount = c(1:12)

# DESeq2 analysis on my local machine
colData = data.frame(condition = factor(c("h48","h48","Input","h0","Input",
                                          "h48","Input","h0","h0")))

# Import the counts into R
all_amp_data = NULL
for (amp in ampcount) {
  temp = read_tsv(paste("Counts/Counts_amp",amp,".txt", sep = '')) %>% 
    select(!"...1") %>% mutate(mutations = Variants) %>%
    separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
    separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
    mutate(Pos = as.double(Pos))
  all_amp_data = rbind(all_amp_data,temp)
}

# do not Remove the nonsense variants
all_amp_data_processed = all_amp_data %>%# filter(Mut != c("X","B")) %>%
  column_to_rownames("Variants")

# Run DeSeq2
countData = all_amp_data_processed %>% select(!c("AA","Pos","Mut")) %>% 
  as.matrix()

dds=DESeqDataSetFromMatrix(countData, colData, formula(~ condition))
dds$condition <- relevel(dds$condition, ref = "Input")
dds=DESeq(dds)

vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)

pca <- prcomp(t(assay(vsd)[,]))
row.names(pca$x)<-make.unique(row.names(pca$x))

condition_names = resultsNames(dds)


# Write out DESeq results against total data 
for (condition in condition_names) {
  res = results(dds, name = condition) %>% 
    as.data.frame() %>%
    rownames_to_column() %>% rename(Variant = rowname)
  write_tsv(res, paste(condition,"_total.txt", sep=""))
}
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup=c("condition"))+theme_classic() #This is the PCA of the aggregate data against all amplicons

ggsave("PCA_plot.png", device = "png")

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Ensure 'condition' is a factor with levels in the desired order
pcaData$condition <- factor(pcaData$condition, levels = c("Input", "h0", "h48"))

write_tsv(pcaData, "pcaData.xls")


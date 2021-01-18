#########################################
# Differential Gene Expression Analysis #
#########################################

## Loading required libraries
library(DESeq2)
library(EnhancedVolcano)
library(gprofiler2)
library(tidyverse)
library(VennDiagram)

## Data import
countdata <- read.csv("data/geneCounts.csv")
countdata <- countdata[, c(2, 7:15)]
colnames(countdata)[1] <- "Geneid"

exclude_genes <- read.csv("data/exclude_genes.csv")

## Remove rRNA/tRNA/noncoding/pseudogenes
countdata <- countdata %>%
  anti_join(exclude_genes, by = "Geneid")
rownames(countdata) <- countdata[,1]
countdata[,1] <- NULL

## Data preprocessing
countdata <- countdata %>%
  select(fex10_1:S2R_3) %>%
  as.matrix()

## Assign condition
(genotype <- factor(c(rep("fex10", 3), rep("fex2", 3), rep("control", 3))))

## Coldata dataframe
(coldata <- data.frame(row.names=colnames(countdata), genotype))

## Instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countdata, 
                              colData=coldata, 
                              design=~genotype)

## Pre-filtering based on read counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

## Retrieve results

# 2.5% Fly Extract vs Control
res_1 <- results(dds, contrast = c("genotype", "fex2", "control"))
resdata_1 <- merge(as.data.frame(res_1), 
                 as.data.frame(counts(dds, normalized=TRUE)), 
                 by="row.names", 
                 sort=FALSE)
names(resdata_1)[1] <- "Symbol"

# 10% Fly Extract vs Control
res_2 <- results(dds, contrast = c("genotype", "fex10", "control"))
resdata_2 <- merge(as.data.frame(res_2), 
                   as.data.frame(counts(dds, normalized=TRUE)), 
                   by="row.names", 
                   sort=FALSE)
names(resdata_2)[1] <- "Symbol"

## Results arrangement and export
resdata_1 <- resdata_1 %>%
  as_tibble() %>%
  arrange(padj) %>% write_csv(file = "data/DESeq2_Results_fex2.csv")

resdata_2 <- resdata_2 %>%
  as_tibble() %>%
  arrange(padj) %>% write_csv(file = "data/DESeq2_Results_fex10.csv")

## Export significantly upregulated and downregulated genes
up_fex2 <- resdata_1 %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Upregulated_genes_fex2.csv")

down_fex2 <- resdata_1 %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Downregulated_genes_fex2.csv")

up_fex10 <- resdata_2 %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Upregulated_genes_fex10.csv")

down_fex10 <- resdata_2 %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Downregulated_genes_fex10.csv")

###########################
# PCA Plot viasualization #
###########################

## Transform count data using variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(dds)

## Generate first two PCs, PC1 and PC2
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## PCA plot visualization with ggplot2
pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position="right", legend.text = element_text(size = 14), 
        legend.title = element_text(face = "bold", size = 16))

pdf("data/PCA_plot.pdf", width=7, height = 6)
pcaplot
dev.off()

#########################
# Intersample variances #
#########################

## Transform count data using regularized logarithm (rlog)
expmatrix_DESeq <- rlog(dds, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)

## Control correlogram and export
ctrlCorr <- expmatrix[,grep("S2R", colnames(expmatrix))]
pdf("data/Ctrl_Correlogram_S2R.pdf", width=5, height=5)
corrgram::corrgram(ctrlCorr, order=TRUE, lower.panel=corrgram::panel.pie,
                   upper.panel=corrgram::panel.pts, text.panel=corrgram::panel.txt,
                   main="Correlogram of Controls",
                   col.regions=colorRampPalette("#1E90FF"))
dev.off()

## Experimental correlogram and export
expCorr_fex2 <- expmatrix[,grep("fex2", colnames(expmatrix))]
pdf("data/Exp_Correlogram_fex2.pdf", width=5, height=5)
corrgram::corrgram(expCorr_fex2, order=TRUE, lower.panel=corrgram::panel.pie,
                   upper.panel=corrgram::panel.pts, text.panel=corrgram::panel.txt,
                   main="Correlogram of 2.5% Fly Extract",
                   col.regions=colorRampPalette("#DC143C"))
dev.off()

expCorr_fex10 <- expmatrix[,grep("fex10", colnames(expmatrix))]
pdf("data/Exp_Correlogram_fex10.pdf", width=5, height=5)
corrgram::corrgram(expCorr_fex10, order=TRUE, lower.panel=corrgram::panel.pie,
                   upper.panel=corrgram::panel.pts, text.panel=corrgram::panel.txt,
                   main="Correlogram of 10% Fly Extract",
                   col.regions=colorRampPalette("#DAA520"))
dev.off()

##########################################
# Gene Ontology Analysis with gProfiler2 #
##########################################

go_up_fex2 <- gost(up_fex2$Symbol, organism = "dmelanogaster", evcodes = TRUE,
              correction_method = "fdr") %>% .$result %>%
  select(c(source, term_id, term_name, p_value, 
           intersection_size, intersection)) %>%
  filter(p_value < 0.05) %>%
  write_csv("data/go_upregulated_fex2.csv")

go_down_fex2 <- gost(down_fex2$Symbol, organism = "dmelanogaster", evcodes = TRUE,
                correction_method = "fdr") %>% .$result %>%
  select(c(source, term_id, term_name, p_value, 
           intersection_size, intersection)) %>%
  filter(p_value < 0.05) %>%
  write_csv("data/go_downregulated_fex2.csv")

##############################
# Volcano plot visualization #
##############################

## Enhanced Volcano Function
EVolcano_Plotter <- function(input_data, p_cutOFF, fc_cutOFF, cus_labels, plt_title){
  EnhancedVolcano(input_data,
                  lab = as.character(input_data$Symbol),
                  x = "log2FoldChange",
                  y = "padj",
                  xlab = bquote(~Log[2]~ "Fold Change"),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  pCutoff = p_cutOFF,
                  FCcutoff = fc_cutOFF,
                  selectLab = cus_labels,
                  drawConnectors = TRUE,
                  colConnectors = "grey30",
                  labCol = "black",
                  labSize = 6,
                  colAlpha = 1,
                  legendPosition = "None",
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  pointSize = 2,
                  title = plt_title,
                  subtitle = "Volcano Plot",
                  caption = paste("adjP cutoff = ", p_cutOFF, " and ", "FC CutOff = ", fc_cutOFF),
                  border = "full",
                  col=c("black", "darkgreen", "goldenrod2", "brown3")
  )
}

## Calling Plotter Function
# Volcano Plot fex2 vs control
cus_labels_1 <- c("PGRP-SA", "CecA2", "Dro", "DptB", "Drs", "upd3", "Pvf3", "Pvf2", "PPO1", "lz",
                  "trol", "Ser", "Glt", "Pxn", "peb", "Tig", "atilla", "Itgbn")
plot_title_1 <- "2.5% Fly Extract vs Control"
vol_1 <- EVolcano_Plotter(resdata_1, 0.05, 1.0, cus_labels_1, plot_title_1)

# Volcano Plot fex10 vs control
cus_labels_2 <- c("PGRP-SA", "CecA2", "Dro", "DptB", "Drs", "upd3", "Pvf3", "Pvf2", "PPO1", "lz", "Egfr", "hh",
                  "trol", "Ser", "Glt", "Pxn", "peb", "Tig", "atilla", "Itgbn")
plot_title_2 <- "10% Fly Extract vs Control"
vol_2 <- EVolcano_Plotter(resdata_2, 0.05, 1.0, cus_labels_2, plot_title_2)

## Saving plots
pdf("data/Volcano_plot_fex2.pdf", width=7, height = 9)
vol_1
dev.off()

pdf("data/Volcano_plot_fex10.pdf", width=7, height = 9)
vol_2
dev.off()
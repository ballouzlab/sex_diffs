#!/usr/bin/env Rscript

# Library path
.libPaths("/directflow/SCCGGroupShare/projects/SeyhanYazar/R_libs_3.6.1")
library("dplyr")
library("valr")
library("data.table")
library("tidyverse")
library("broom")
library("ggplot2")
# library("vroom")

# Call the cell types
args = commandArgs(trailingOnly=TRUE)

genename <- args[2]
rsid <- args[1]

# rsid <- "rs231770"
# genename <- "CTLA4"

# Colour palettee
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", 
               "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
# Order of independent eQTLs
rank_order <- c("eSNP1","eSNP2","eSNP3", "eSNP4","eSNP5")

# Order of cells
cell_order <- c("CD4all","CD4effCM","CD4TGFbStim","CD8eff","CD8all","CD8unknown",
                "NKmat","NKact","Plasma","Bmem","BimmNaive","MonoC","MonoNC","DC")

cell_labels <- c("CD4 NC","CD4 ET","CD4 SOX4","CD8 ET","CD8 NC","CD8 S100B","NK", "NK R",
                 "Plasma","B Mem","B IN","Mono C","Mono NC", "DC")

cells <- data.frame(cell_order,cell_labels)
colnames(cells)[1] <- "cellType"
cells$cellType <- as.character(cells$cellType)

df <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/main_analysis/genotype_expression_plots/%s_%s_residuals.tsv", genename, rsid))
df <- left_join(df, cells, by="cellType")
colnames(df)[2] <- "residuals"
df$cell_labels <- factor(df$cell_labels, levels=cell_labels)

gp <- ggplot(df, aes(x = genotype, y = residuals, color=cell_labels)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ cell_labels, scales = "free", nrow=2) + 
  theme_classic() +
  scale_color_manual(values=tol14rainbow, drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=cell_labels), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))

ggsave(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/main_analysis/genotype_expression_plots/%s_%s_residual_plot_211011.pdf", genename, rsid),  
    width=10, height=4.0)
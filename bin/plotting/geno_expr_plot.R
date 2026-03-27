library("dplyr")
library("valr")
library("data.table")
library("tidyverse")
library("broom")
library("ggplot2")

setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-31_sex-specific-eQTL-expr-geno-plots")

male <- fread("TCF4_rs9952839_plot_data.tsv") %>% filter(cellType=="B_memory")
genename <- "TCF4"
rsid <- "rs9952839"
male$sex <- factor(male$sex, levels=c("joint", "female", "male"))

color_palette <- c("#2A788EFF", "#671D5B","#F7A51D")

gp <- ggplot(male, aes(x = genotype, y = expression, color=sex)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ sex, scales = "free", nrow=1) + 
  theme_classic(18) +
  scale_color_manual(values=color_palette, drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=sex), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))

ggsave(sprintf("%s_%s_residual_plot_211011.pdf", genename, rsid),  
    width=8.0, height=4.0)


male <- fread("ICAM3_rs8104608_plot_data.tsv") %>% filter(cellType=="B_memory")
genename <- "ICAM3"
rsid <- "rs8104608"
male$sex <- factor(male$sex, levels=c("joint", "female", "male"))

color_palette <- c("#2A788EFF", "#671D5B","#F7A51D")

gp <- ggplot(male, aes(x = genotype, y = expression, color=sex)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ sex, scales = "free", nrow=1) + 
  theme_classic(18) +
  scale_color_manual(values=color_palette, drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=sex), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))

ggsave(sprintf("%s_%s_residual_plot.pdf", genename, rsid),  
    width=8.0, height=4.0)

male <- fread("ACP5_rs117694210_plot_data.tsv") %>% filter(cellType=="B_memory")
genename <- "ACP5"
rsid <- "rs117694210"
male$sex <- factor(male$sex, levels=c("joint", "female", "male"))

color_palette <- c("#2A788EFF", "#671D5B","#F7A51D")

gp <- ggplot(male, aes(x = genotype, y = expression, color=sex)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ sex, scales = "free", nrow=1) + 
  theme_classic(18) +
  scale_color_manual(values=color_palette, drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=sex), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))

ggsave(sprintf("%s_%s_residual_plot.pdf", genename, rsid),  
    width=8.0, height=4.0)

male <- fread("COX20_rs9730334_plot_data.tsv ") %>% filter(cellType=="B_memory")
genename <- "COX20"
rsid <- "rs9730334"
male$sex <- factor(male$sex, levels=c("joint", "female", "male"))

color_palette <- c("#2A788EFF", "#671D5B","#F7A51D")

gp <- ggplot(male, aes(x = genotype, y = expression, color=sex)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ sex, scales = "free", nrow=1) + 
  theme_classic(18) +
  scale_color_manual(values=color_palette, drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=sex), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))

ggsave(sprintf("%s_%s_residual_plot.pdf", genename, rsid),  
    width=8.0, height=4.0)
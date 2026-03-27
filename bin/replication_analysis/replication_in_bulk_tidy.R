##############################################################################
# Script information                                                      
# Title: sex specific eqtl analysis
# Author: Seyhan Yazar
# Date: 2020-10-18
# Description: 
##############################################################################

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-11-24", "bulk-replication")

library(data.table)
library(dplyr)

load("data/overlaps_sex_egenes.Rdata")

df <- temp %>% select("geneid","snpid")
fwrite(df, "data/sex_specific_examples.txt", sep=" ", col.names=FALSE)
df_f <- temp[,c(16,15)]
fwrite(df_f, "data/sex_specific_examples_f.txt", sep=" ", col.names=FALSE)

files <- list.files(path=output, pattern="porcu_replication_*", full.names=T)
rep_list <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
rep_df <- bind_rows(rep_list)
dim(rep_df)

rep_df$significant_joint <- ifelse(rep_df$localFDR_j_onek1k < 0.05,1,0)
rep_df$significant_m <- ifelse(rep_df$localFDR_m_onek1k < 0.05,10,0)
rep_df$significant_f <- ifelse(rep_df$localFDR_f_onek1k < 0.05,100,0)
rep_df$significance <- rep_df$significant_joint + rep_df$significant_m + rep_df$significant_f
rep_df$significance[is.na(rep_df$significance)] <- "not_tested"
rep_df$significance <- factor(rep_df$significance, labels=c("not_significant_in_all",
    "significant_only_in_joint",
    "significant_in_joint_and_male",
    "significant_in_joint_and_female",
    "significant_in_all", 
    "not_tested"))
table(rep_df$significance)

library(ggplot2)
library(viridis)
library(hrbrthemes)

rep_barplot <- rep_df %>% 
    ggplot(.,aes(fill=significance, x=celltype)) + 
    geom_bar(position="fill", stat="count") + 
    scale_fill_viridis(discrete = T, name="") +
    ggtitle("Replication of sex-interacting eQTLs from Porcu et al (n=460)") +
    theme_ipsum() +
    xlab("") + ylab("Percentage(%)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sprintf("%s/porcu_replication_eQTL_barplot.png",output), rep_barplot, 
    width=10, height=7)

fwrite(rep_df, sprintf("%s/porcu_replication_results.tsv",output), sep="\t")

##### Replication Olivia et al 2022

files <- list.files(path=output, pattern="olivia_replication_*", full.names=T)
rep_list <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
rep_df <- bind_rows(rep_list)
dim(rep_df)

rep_df$significant_joint <- ifelse(rep_df$localFDR_j_onek1k < 0.05,1,0)
rep_df$significant_m <- ifelse(rep_df$localFDR_m_onek1k < 0.05,10,0)
rep_df$significant_f <- ifelse(rep_df$localFDR_f_onek1k < 0.05,100,0)
rep_df$significance <- rep_df$significant_joint + rep_df$significant_m + rep_df$significant_f
rep_df$significance[is.na(rep_df$significance)] <- "not_tested"
rep_df$significance <- factor(rep_df$significance, labels=c("not_significant_in_any",
    "significant_only_in_joint",
    "significant_in_joint_and_female",
    "significant_in_all", 
    "not_tested"))
table(rep_df$significance)

library(ggplot2)
library(viridis)
library(hrbrthemes)

rep_barplot <- rep_df %>% 
    ggplot(.,aes(fill=significance, x=celltype)) + 
    geom_bar(position="fill", stat="count") + 
    scale_fill_viridis(discrete = T, name="") +
    ggtitle("Replication of sex-interacting eQTLs from Olivia et al (n=3)") +
    theme_ipsum() +
    xlab("") + ylab("Percentage(%)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(sprintf("%s/olivia_replication_eQTL_barplot.png",output), rep_barplot, 
    width=10, height=7)

fwrite(rep_df, sprintf("%s/olivia_replication_results.tsv",output), sep="\t")

##### Replication Kukarba

files <- list.files(path=output, pattern="kukarba_replication_*", full.names=T)
rep_list <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
rep_df <- bind_rows(rep_list)
dim(rep_df)

rep_df$significant_joint <- ifelse(rep_df$localFDR_j_onek1k < 0.05,1,0)
rep_df$significant_m <- ifelse(rep_df$localFDR_m_onek1k < 0.05,10,0)
rep_df$significant_f <- ifelse(rep_df$localFDR_f_onek1k < 0.05,100,0)
rep_df$significance <- rep_df$significant_joint + rep_df$significant_m + rep_df$significant_f
rep_df$significance[is.na(rep_df$significance)] <- "not_tested"
rep_df$significance <- factor(rep_df$significance, labels=c("not_significant_in_any",
    "not_tested"))
table(rep_df$significance)

fwrite(rep_df, sprintf("%s/kukarba_replication_results.tsv",output), sep="\t")

##### Replication Yao

files <- list.files(path=output, pattern="yao_replication_*", full.names=T)
rep_list <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
rep_df <- bind_rows(rep_list)
dim(rep_df)

rep_df$significant_joint <- ifelse(rep_df$localFDR_j_onek1k < 0.05,1,0)
rep_df$significant_m <- ifelse(rep_df$localFDR_m_onek1k < 0.05,10,0)
rep_df$significant_f <- ifelse(rep_df$localFDR_f_onek1k < 0.05,100,0)
rep_df$significance <- rep_df$significant_joint + rep_df$significant_m + rep_df$significant_f
rep_df$significance[is.na(rep_df$significance)] <- "not_tested"
rep_df$significance <- factor(rep_df$significance, labels=c("not_significant_in_any",
    "not_tested"))
table(rep_df$significance)

fwrite(rep_df, sprintf("%s/yao_replication_results.tsv",output), sep="\t")

### Collate all studies ####

porcu <- fread(sprintf("%s/porcu_replication_results.tsv",output)) %>% select("celltype","SNP", "geneid","snpid",14:23,"significance")
porcu$study <- "Porcu et al 2022"
olivia <- fread(sprintf("%s/olivia_replication_results.tsv",output)) %>% select("celltype","SNP", "geneid","snpid",26:34,"significance")
olivia$study <- "Olivia et al 2020"
kukarba <- fread(sprintf("%s/kukarba_replication_results.tsv",output)) %>% select("celltype","SNP", "geneid","snpid",14:22,"significance")
kukarba$study <- "Kukarba et al 2016"
yao <- fread(sprintf("%s/yao_replication_results.tsv",output)) 
yao_df <- fread("data/bulk_data/Yao_sex_interacting_eQTLs_with_rsID.csv") %>% select("ID","snpid")
yao <- left_join(yao, yao_df, by= "ID")
colnames(yao) [1] <- "SNP"
yao <- yao %>% select("celltype","SNP", "geneid","snpid",10:18,"significance")
yao$study <- "Yao et al 2014"

rep_studies <- bind_rows(porcu, olivia,kukarba, yao)
rep_studies2 <- rep_studies %>% 
    filter(!significance==c("not_tested")) 

new_vars <- list("not_significant_in_all" = 0,
    "not_significant_in_any" = 0,
    "significant_only_in_joint" = 1, 
    "significant_in_joint_and_female" = 1,
    "significant_in_joint_and_male" = 1,
    "significant_in_all" = 1)

rep_studies2 <- rep_studies2 %>% mutate(signif = new_vars[significance])
rep_studies2$signif <- unlist(rep_studies2$signif)

rep_studies2 %>% group_by(study, celltype) %>% count()

library(ggplot2)

rep_studies2 %>% 
    group_by(study, geneid, celltype) %>% 
    summarise(no_snps = sum(signif)) %>%
    mutate(replication=ifelse(no_snps==0, 0,1)) %>% 
    group_by(study, geneid) %>%
    summarise(no_celltypes=sum(replication)) %>%
    ggplot(.,aes(study, geneid, fill = no_celltypes)) +
    geom_tile() +
    theme_bw()

ggsave(file="results/Figures_and_Tables/replication_plots/porcu_test.png")

rep_studies2$study_f <- factor(rep_studies2$study, levels=c("Yao et al 2014","Kukarba et al 2016","Olivia et al 2020", "Porcu et al 2022"))

rep_studies2 %>%
    group_by(study_f, geneid, celltype) %>% 
    summarise(no_snps = sum(signif)) %>%
    mutate(replication=ifelse(no_snps==0, 0,1)) %>%
    mutate_at(vars(replication), factor) %>% 
    ggplot(., aes(x = celltype, y = geneid, fill = replication)) +
    geom_tile(colour="white", size=0.50, stat="identity") +
    scale_fill_manual(values=c("grey30","steelblue"), labels=c("not replicated","replicated at least in one analysis")) +
    facet_wrap(~study_f,nrow=1) + 
    theme_gray() +
    theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
    labs(y="Gene") +
    theme(axis.title.x = element_blank(), legend.position="bottom", legend.title = element_blank()) 

ggsave(file="results/Figures_and_Tables/replication_plots/heatmap_replication.pdf", width=12, height=7)

ggsave(file="results/Figures_and_Tables/replication_plots/heatmap_test.png", width=12, height=7)




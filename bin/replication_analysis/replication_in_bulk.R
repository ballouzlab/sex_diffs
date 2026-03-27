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

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)

# Example
# cellLabel <- "B_memory"

library(data.table)
library(dplyr)

# Read HRC ids and add to the dataframes
# hrc_ids <- fread("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
# head(hrc_ids)
# hrc_ids$snpid <- sprintf("%s:%s_%s", hrc_ids$'#CHROM', hrc_ids$POS, hrc_ids$ALT)
# hrc_ids <- hrc_ids %>% select(snpid, '#CHROM', POS, ID)
# colnames(hrc_ids) <- c("snpid","chr","pos","SNP")

# porcu_df <- fread("data/bulk_data/Porcu_sex_biased_eQTLs.csv")
# porcu_df <- left_join(porcu_df, hrc_ids, by="SNP")
# table(porcu_df$chr)
# fwrite(porcu_df,"data/bulk_data/Porcu_sex_biased_eQTLs_with_rsID.csv")

# olivia_df <- fread("data/bulk_data/Olivia_2020_sex_biased_eQTLs.tsv")
# colnames(olivia_df)[5] <- "SNP"
# olivia_df <- left_join(olivia_df,hrc_ids, by="SNP")
# fwrite(olivia_df,"data/bulk_data/Olivia_sex_biased_eQTLs_with_rsID.csv")

# porcu_df <- fread("data/bulk_data/Porcu_sex_biased_eQTLs_with_rsID.csv")
# 
# chr_to_look <- c(unique(porcu_df$chr)[!is.na(unique(porcu_df$chr))])
# 
# extract_joint_results <- function(x) {
#     chr_number <- x
#     bulk <- porcu_df %>% filter(chr==chr_number)
#     colnames(bulk)[3] <- "geneid"
#     joint <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/%s/round1/%s_chr%s_correlation_results.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, joint, by=c("snpid", "geneid"))
#     male <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr%s_correlation_results_male.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, male, by=c("snpid", "geneid"))
#     female <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/%s/round1/%s_chr%s_correlation_results_female.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, female, by=c("snpid", "geneid"))
#     bulk
# }
# 
# replication_results <- lapply(chr_to_look, extract_joint_results)
# replication_results_df <- bind_rows(replication_results)
# colnames(replication_results_df)
# replication_results_df <- replication_results_df %>% select(c("SNP","Gene","geneid","Allele",
#     "zscore_f", "zscore_m", "Z_diff", "pval",
#     "FDR", "pvalue_f", "pvalue_m", "snpid",
#     "chr", "pos", "estimate.x", "p.value.x","localFDR.x",
#     "estimate.y", "p.value.y", "localFDR.y",
#     "estimate", "p.value", "localFDR"))
# colnames(replication_results_df)
# colnames(replication_results_df)[15:23] <- c("estimate_j_onek1k", "pvalue_j_onek1k", "localFDR_j_onek1k",
#     "estimate_m_onek1k", "pvalue_m_onek1k", "localFDR_m_onek1k",
#     "estimate_f_onek1k", "pvalue_f_onek1k", "localFDR_f_onek1k")
# replication_results_df$celltype <- cellLabel
# head(replication_results_df)
# 
# fwrite(replication_results_df, sprintf("%s/porcu_replication_%s.tsv", output, cellLabel),sep="\t")

##### Replication for Olivia et al 2020 #####

# olivia_df <- fread("data/bulk_data/Olivia_sex_biased_eQTLs_with_rsID.csv")
# colnames(olivia_df)[2] <- "geneid"
# olivia_df <- olivia_df %>% filter(Tissue=="Whole_Blood")
# 
# chr_to_look <- c(unique(olivia_df$chr)[!is.na(unique(olivia_df$chr))])
# extract_joint_results <- function(x) {
#     chr_number <- x
#     bulk <- olivia_df %>% unique() %>% filter(chr==chr_number)
#     # colnames(bulk)[3] <- "geneid"
#     joint <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/%s/round1/%s_chr%s_correlation_results.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, joint, by=c("snpid", "geneid"))
#     male <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr%s_correlation_results_male.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, male, by=c("snpid", "geneid"))
#     female <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/%s/round1/%s_chr%s_correlation_results_female.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, female, by=c("snpid", "geneid"))
#     bulk
# }
# 
# replication_results <- lapply(chr_to_look[!chr_to_look %in% c("X","")], extract_joint_results)
# 
# replication_results_df <- bind_rows(replication_results)
# colnames(replication_results_df)
# replication_results_df <- replication_results_df %>% select(c(1:25, "estimate.x", "p.value.x","localFDR.x",
#     "estimate.y", "p.value.y", "localFDR.y",
#     "estimate", "p.value", "localFDR"))
# colnames(replication_results_df)
# colnames(replication_results_df)[26:34] <- c("estimate_j_onek1k", "pvalue_j_onek1k", "localFDR_j_onek1k",
#     "estimate_m_onek1k", "pvalue_m_onek1k", "localFDR_m_onek1k",
#     "estimate_f_onek1k", "pvalue_f_onek1k", "localFDR_f_onek1k")
# replication_results_df$celltype <- cellLabel
# head(replication_results_df)
# 
# fwrite(replication_results_df, sprintf("%s/olivia_replication_%s.tsv", output, cellLabel),sep="\t")

##### Replication for Kukarba 2016 #####
# kukarba_df <- fread("data/bulk_data/Kukarba_sex_interacting_eQTLs.tsv")
# colnames(kukarba_df)[4] <- "geneid"
# 
# chr_to_look <- c("5","16","3")
# extract_joint_results <- function(x) {
#     chr_number <- x
#     bulk <- kukarba_df %>% filter(CHR==chr_number)
#     # colnames(bulk)[3] <- "geneid"
#     joint <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/%s/round1/%s_chr%s_correlation_results.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, joint, by=c("snpid", "geneid"))
#     male <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr%s_correlation_results_male.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, male, by=c("snpid", "geneid"))
#     female <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/%s/round1/%s_chr%s_correlation_results_female.tsv", cellLabel, cellLabel,chr_number))
#     bulk <- left_join(bulk, female, by=c("snpid", "geneid"))
#     bulk
# }
# 
# replication_results <- lapply(chr_to_look, extract_joint_results)
# replication_results_df <- bind_rows(replication_results)
# colnames(replication_results_df)
# replication_results_df <- replication_results_df %>% select(c(1:13, "estimate.x", "p.value.x","localFDR.x",
#     "estimate.y", "p.value.y", "localFDR.y",
#     "estimate", "p.value", "localFDR"))
# colnames(replication_results_df)
# colnames(replication_results_df)[14:22] <- c("estimate_j_onek1k", "pvalue_j_onek1k", "localFDR_j_onek1k",
#     "estimate_m_onek1k", "pvalue_m_onek1k", "localFDR_m_onek1k",
#     "estimate_f_onek1k", "pvalue_f_onek1k", "localFDR_f_onek1k")
# replication_results_df$celltype <- cellLabel
# head(replication_results_df)
# 
# fwrite(replication_results_df, sprintf("%s/kukarba_replication_%s.tsv", output, cellLabel),sep="\t")

##### Replication for Yao 2014 #####
# yao_df <- fread("data/bulk_data/Yao_2014_sex_interacting_eQTLs.csv")
# colnames(yao_df)[1] <- "ID"
# yao_df <- left_join(yao_df, hrc_ids, by="ID")
# yao_df$snpid <- sprintf("%s:%s_%s", yao_df$'#CHROM', yao_df$POS, yao_df$ALT)
# colnames(yao_df)[10] <- "CHR"
# fwrite(yao_df,"data/bulk_data/Yao_sex_interacting_eQTLs_with_rsID.csv")

yao_df <- fread("data/bulk_data/Yao_sex_interacting_eQTLs_with_rsID.csv")

chr_to_look <- unique(yao_df$CHR)
extract_joint_results <- function(x) {
    chr_number <- x
    bulk <- yao_df %>% filter(CHR==chr_number)
    # colnames(bulk)[3] <- "geneid"
    joint <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/results/2021-08-16_conditional-analysis/%s/round1/%s_chr%s_correlation_results.tsv", cellLabel, cellLabel,chr_number))
    bulk <- left_join(bulk, joint, by=c("snpid", "geneid"))
    male <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr%s_correlation_results_male.tsv", cellLabel, cellLabel,chr_number))
    bulk <- left_join(bulk, male, by=c("snpid", "geneid"))
    female <- fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/%s/round1/%s_chr%s_correlation_results_female.tsv", cellLabel, cellLabel,chr_number))
    bulk <- left_join(bulk, female, by=c("snpid", "geneid"))
    bulk
}

replication_results <- lapply(chr_to_look, extract_joint_results)
replication_results_df <- bind_rows(replication_results)
colnames(replication_results_df)
replication_results_df <- replication_results_df %>% select(c(1:9, "estimate.x", "p.value.x","localFDR.x",
    "estimate.y", "p.value.y", "localFDR.y",
    "estimate", "p.value", "localFDR"))
colnames(replication_results_df)
colnames(replication_results_df)[10:18] <- c("estimate_j_onek1k", "pvalue_j_onek1k", "localFDR_j_onek1k",
    "estimate_m_onek1k", "pvalue_m_onek1k", "localFDR_m_onek1k",
    "estimate_f_onek1k", "pvalue_f_onek1k", "localFDR_f_onek1k")
replication_results_df$celltype <- cellLabel
head(replication_results_df)

fwrite(replication_results_df, sprintf("%s/yao_replication_%s.tsv", output, cellLabel),sep="\t")
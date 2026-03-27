#!/usr/bin/env Rscript

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-10-18", "pi-analysis")

# Upload libraries
library(data.table)
library(dplyr)
library(stringr)
library(qvalue)

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
celltype1 <- args[1]

# celltype1 <- "B_intermediate"
out.dir = paste0(output,'/p_outputs/')
dir.create(out.dir)

# Get all the result files (significant eQTls)
female_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/female_eqtls/female_top_eqtls.tsv"
male_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/male_eqtls/male_top_eqtls.tsv"
joint_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_joint_analysis_one_round/joint_analysis_one_round_results.tsv"

df_female <- fread(female_file)
df_male <- fread(male_file)
df_joint <- fread(joint_file)
df_joint <- df_joint %>% 
    mutate(snpid=paste0(Chromosome,":",Position,"_", SNP_assessed_allele))
colnames(df_joint) [1:2] <- c("celltype", "geneid")

female_results_df <- df_female %>% filter(celltype==celltype1)
male_results_df <- df_male %>% filter(celltype==celltype1)
joint_results_df <- df_joint %>% filter(celltype==celltype1)

# Exract female results that are not in the joint then male outputs
not_in_joint_df <- female_results_df[!female_results_df$snpid %in% joint_results_df$snpid,]
dim(not_in_joint_df)
not_in_male_df <- not_in_joint_df[!not_in_joint_df$snpid %in% male_results_df$snpid,]
not_in_male_df <- not_in_male_df %>% select(geneid, snpid)
dim(not_in_male_df)
head(not_in_male_df)

extract_info <- function (x,y) {
    GENE=y$geneid
    SNP=y$snpid
    # GENE="PPA1"
    # SNP="10:71993094_G"
    chr=sub(":.*","",SNP)
    filename <- sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-09-27_male-only/%s/round1/%s_chr%s_correlation_results_male.tsv", celltype1, celltype1, chr)
    df_original_eqtl <- fread(filename)
    df_original_eqtl <- df_original_eqtl %>% filter(snpid==SNP & geneid==GENE)
    df_original_eqtl %>% select(-c("geneid", "snpid"))
}

pvalues_df <- not_in_male_df %>% group_by(snpid,geneid) %>% group_modify(extract_info)

fwrite(pvalues_df, sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/p_outputs/female_vs_male_%s_pvalues.tsv", celltype1), quote=F, sep="\t")

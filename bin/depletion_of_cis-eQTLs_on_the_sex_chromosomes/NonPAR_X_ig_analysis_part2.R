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
output <- set_output("2022-12-09", "sex-chr-analysis")

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)

# Example
# cellLabel <- "dnT"

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)

input.dir = paste0(output,'/NonPAR-X-ig/', cellLabel)

# Input files
log_residual_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_X_ig_log_residuals.tsv'))
genotype_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_X_ig_genotypes.tsv'))
gene_snp_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_X_ig_gene-snp_pairs.tsv'))

# mann_whitney's u-test to compare independent groups
# x is the data frame with chr and pos, y is snpid and geneid
mann_whitney <- function (x,y) {
  gene <- y$geneid
  # print(gene)
  snp <- y$snpid
  # print(snp)

  # Select values to test
  res_val <- log_residual_df %>% select(c("sampleid",all_of(gene)))
  # class(res_val)
  genotype_val <- genotype_df[,c("sampleid",snp), with = FALSE] 

  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by=c("sampleid"))
  colnames(test_df) <- c("sampleid","residual", "SNP")
  
  # Generate model
  # model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=T)
  # model_table <- tidy(model)
  # model_table
  
  model2 <- wilcox.test(residual~SNP, data=test_df)
  model_table2 <- tidy(model2)
  model_table2
}

gene_snp_df <- gene_snp_df %>% filter(snpid %in% colnames(genotype_df))

gene_snp_df_splitted <- split(gene_snp_df, sample(1:5, nrow(gene_snp_df), replace=T))

for (i in 1:3) {
    gene_snp_df_part <- gene_snp_df_splitted[i] %>% data.frame()
    colnames(gene_snp_df_part) <-  c("chr","pos","snpid","geneid")
    log_spearman_df <- gene_snp_df_part %>% group_by(geneid,snpid) %>% group_modify(mann_whitney)
    fwrite(log_spearman_df,sprintf("%s/%s_nonPAR_X_ig_correlation_results_%s.tsv", input.dir, cellLabel,i),sep="\t",quote=F)
}
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

input.dir = paste0(output,'/NonPAR-XaXiX/', cellLabel)

# Input files
log_residual_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_XaXiX_log_residuals.tsv'))
genotype_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_XaXiX_genotypes.tsv'))
gene_snp_df <- fread(paste0(input.dir,"/", cellLabel,'_nonPAR_XaXiX_gene-snp_pairs.tsv'))

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
spearman_correlation <- function (x,y) {
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
  model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=T)
  model_table <- tidy(model)
  model_table
}

gene_snp_df_splitted <- split(gene_snp_df, sample(1:5, nrow(gene_snp_df), replace=T))

for (i in 1:5) {
    gene_snp_df_part <- gene_snp_df_splitted[i] %>% data.frame()
    colnames(gene_snp_df_part) <-  c("chr","pos","snpid","geneid")
    log_spearman_df <- gene_snp_df_part %>% group_by(geneid,snpid) %>% group_modify(spearman_correlation)
    fwrite(log_spearman_df,sprintf("%s/%s_nonPAR_XaXiX_correlation_results_%s.tsv", input.dir, cellLabel,i),sep="\t",quote=F)
}

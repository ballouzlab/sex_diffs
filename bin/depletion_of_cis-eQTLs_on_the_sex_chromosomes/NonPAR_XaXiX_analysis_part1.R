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
# cellLabel <- "B_memory"
# cellLabel <- "dnT"

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)

# Input files
genotype_df <- fread("data/X/chrX_genotypes.tsv")
snpLoc_df <- fread("data/X/chrX.filtered.bim")
geneLoc_df <- fread("data/geneloc_chrX.tsv")
colnames(geneLoc_df)[1] <- "geneid"
escape_genes_df <- fread("data/NonPAR_escape_genes.tsv")
not_escape_genes_df <- geneLoc_df %>% filter(!geneid %in% escape_genes_df$geneid)
geneLoc_df_subPARX <- fread("data/subPARX.tsv")
not_escape_genes_df <- not_escape_genes_df %>% filter(!geneid %in% geneLoc_df_subPARX$geneid)

data.dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21'
covariate_df <- fread(sprintf("%s/outputs/peer_factors/%s_peer_factors_after_removing_genes.tsv", data.dir, cellLabel))
expression_df <- fread(sprintf("%s/outputs/filtered_matrices/%s_expression_genes_removed.tsv", data.dir, cellLabel))

## Count matrix
dim(expression_df)
expression_df <- expression_df[,-c("V1")]
print(expression_df[1:5,1:5])
sample_ids <- expression_df$sampleid

# Individuals with non-zero expression has been eliminated previously
# log+1 transformation
logplusone <- function(x) {log(x[1] + 1)}
log_exprs_mat <- apply(expression_df[,2:(length(expression_df))],1:2,logplusone)
log_exprs_df <- setDT(as.data.frame(log_exprs_mat))
log_exprs_df$sampleid <- sample_ids
print(log_exprs_df[1:5,1:5])

out.dir = paste0(output,'/NonPAR-XaXiX/')
dir.create(out.dir)
out.dir2 = paste0(out.dir,"/",cellLabel)
dir.create(out.dir2)

# Extract PAR X region genes 
log_exprs_df_parx <- log_exprs_df %>% select("sampleid",colnames(log_exprs_df)[(colnames(log_exprs_df) %in% not_escape_genes_df$geneid) == TRUE])
dim(log_exprs_df_parx)

fwrite(log_exprs_df_parx,paste0(out.dir2,"/", cellLabel,'_nonPAR_XaXiX_expression.tsv'),sep="\t")

dim(snpLoc_df)
colnames(snpLoc_df) <- c('chr', 'snpid','ignore','pos','A1','A2')
head(snpLoc_df)

# Find Gene-SNP pairs for chrNumber
gene_ids <- colnames(log_exprs_df_parx)[-1]
geneLoc_df <- not_escape_genes_df[(not_escape_genes_df$geneid %in% gene_ids),]
geneLoc_df$left <- geneLoc_df$start - 1000000
geneLoc_df$right <- geneLoc_df$end + 1000000
gene_snp_df <- geneLoc_df[!duplicated(geneLoc_df$geneid)] %>% 
  group_by(geneid) %>% 
  group_modify(function(x, y) tibble(snpid = snpLoc_df$snpid[between(as.numeric(snpLoc_df$pos), as.numeric(x["left"]), as.numeric(x["right"]))]), .keep=FALSE)
# gene_snp_df <- gene_snp_df %>% ungroup()
gene_snp_df <- inner_join(snpLoc_df, gene_snp_df) %>% select(chr, pos, snpid, geneid)
dim(gene_snp_df)
print(gene_snp_df %>% head)

fwrite(gene_snp_df,paste0(out.dir2,"/", cellLabel,'_nonPAR_XaXiX_gene-snp_pairs.tsv'),sep="\t")

# Find residuals for log transformed expressions
calculate_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- log_exprs_df_parx %>% select("sampleid",all_of(gene))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,covariate_df, by="sampleid")
  colnames(test_df)[2] <- "expression"

  # Generate model
  model <- lm(expression ~ age + sex + pc1 + pc2 + pc3 + pc4 + pf1 + pf2 , data=test_df)
  residuals=resid(model)
  residuals  
}

log_residual_mat <- sapply(gene_ids,calculate_residuals)
rownames(log_residual_mat) <- log_exprs_df_parx$sampleid
log_residual_df <- data.frame(log_residual_mat, check.names=FALSE)
log_residual_df$sampleid <- log_exprs_df_parx$sampleid
print(log_residual_df[1:5,1:5])
dim(log_residual_df)

fwrite(log_residual_df,paste0(out.dir2,"/", cellLabel,'_nonPAR_XaXiX_log_residuals.tsv'),sep="\t")

genotype_df <- genotype_df %>% filter(sampleid %in% log_exprs_df_parx$sampleid)

fwrite(genotype_df,paste0(out.dir2,"/", cellLabel,'_nonPAR_XaXiX_genotypes.tsv'),sep="\t")

print('JOB IS DONE!')

quit()
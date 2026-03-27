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
# cellLabel <- "MAIT"

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
escape_genes_df <- left_join(escape_genes_df, geneLoc_df, by="geneid")
geneLoc_df_subPARX <- fread("data/subPARX.tsv")
escape_genes_df <- escape_genes_df %>% filter(!geneid %in% geneLoc_df_subPARX$geneid)

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

out.dir = paste0(output,'/NonPAR-X-eg/')
dir.create(out.dir)
out.dir2 = paste0(out.dir,"/",cellLabel)
dir.create(out.dir2)

# Separate covariate files for females

no_males = covariate_df %>% filter(sex==1) %>% nrow() 
no_females = covariate_df %>% filter(sex==2) %>% nrow()

covariate_df_m <- covariate_df %>% filter(sex==1) 
dim(covariate_df_m)

fwrite(covariate_df_m,paste0(out.dir2,'/',cellLabel,'_male_covariates.tsv'),sep="\t")

# Separate expression file for females

log_exprs_df <- expression_df %>% filter(sampleid %in% covariate_df_m$sampleid)
dim(log_exprs_df)

# Extract PAR X region genes 
log_exprs_df_parx <- log_exprs_df %>% select("sampleid",colnames(log_exprs_df)[(colnames(log_exprs_df) %in% escape_genes_df$geneid) == TRUE])
dim(log_exprs_df_parx)

fwrite(log_exprs_df_parx,paste0(out.dir2,"/", cellLabel,'_nonPAR_X_eg_expression.tsv'),sep="\t")

## SNP location file

dim(snpLoc_df)
colnames(snpLoc_df) <- c('chr', 'snpid','ignore','pos','A1','A2')
head(snpLoc_df)

# Find Gene-SNP pairs for chrNumber
gene_ids <- colnames(log_exprs_df_parx)[-1]
geneLoc_df <- escape_genes_df[(escape_genes_df$geneid %in% gene_ids),]
geneLoc_df$left <- geneLoc_df$start - 1000000
geneLoc_df$right <- geneLoc_df$end + 1000000
gene_snp_df <- geneLoc_df[!duplicated(geneLoc_df$geneid)] %>% 
  group_by(geneid) %>% 
  group_modify(function(x, y) tibble(snpid = snpLoc_df$snpid[between(as.numeric(snpLoc_df$pos), as.numeric(x["left"]), as.numeric(x["right"]))]), .keep=FALSE)
# gene_snp_df <- gene_snp_df %>% ungroup()
gene_snp_df <- inner_join(snpLoc_df, gene_snp_df) %>% select(chr, pos, snpid, geneid)
dim(gene_snp_df)
print(gene_snp_df %>% head)

fwrite(gene_snp_df,paste0(out.dir2,"/", cellLabel,'_nonPAR_X_eg_gene-snp_pairs.tsv'),sep="\t")

# Find residuals for log transformed expressions
calculate_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- log_exprs_df_parx %>% select("sampleid",all_of(gene))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,covariate_df_m,by="sampleid")
  colnames(test_df)[2] <- "expression"

  # Generate model
  model <- lm(expression ~ age + pc1 + pc2 + pc3 + pc4 + pf1 + pf2 , data=test_df)
  residuals=resid(model)
  residuals  
}

log_residual_mat <- sapply(gene_ids,calculate_residuals)
rownames(log_residual_mat) <- log_exprs_df_parx$sampleid
log_residual_df <- data.frame(log_residual_mat, check.names=FALSE)
log_residual_df$sampleid <- log_exprs_df_parx$sampleid
print(log_residual_df[1:5,1:5])
dim(log_residual_df)

fwrite(gene_snp_df,paste0(out.dir2,"/", cellLabel,'_nonPAR_X_eg_log_residuals.tsv'),sep="\t")

genotype_df <- genotype_df %>% filter(sampleid %in% log_exprs_df_parx$sampleid)

# Filter any row that has a 1 in the column
# genotype_df %>% filter(across(everything(), ~ !grepl("1", .)))

# Filter any column that has 1
genotype_df %>% select(where(~is.numeric(.x) && any(.x == 1))) %>% ncol()
genotype_df2 <- genotype_df %>% select("sampleid",!where(~is.numeric(.x) && any(.x == 1)))

fwrite(genotype_df2,paste0(out.dir2,"/", cellLabel,'_nonPAR_X_eg_genotypes.tsv'),sep="\t")
# fwrite(genotype_df2,paste0(out.dir2,"/", cellLabel,'_nonPAR_X_eg_genotypes_errors.tsv'),sep="\t")

# mann_whitney's u-test to compare independent groups
# x is the data frame with chr and pos, y is snpid and geneid
mann_whitney <- function (x,y) {
  gene <- y$geneid
  # print(gene)
  snp <- y$snpid
  # print(snp)

  # Select values to test
  res_val <- log_residual_df[c("sampleid",gene)]
  # class(res_val)
  genotype_val <- genotype_df2[,c("sampleid",snp), with = FALSE] 

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

gene_snp_df <- gene_snp_df %>% filter(snpid %in% colnames(genotype_df2))

log_mann_whitney_df <- gene_snp_df %>% group_by(geneid,snpid) %>% group_modify(mann_whitney)

# Calculate the qvalues for pvalues
pvalues <- log_mann_whitney_df$p.value
qobj <- qvalue(p = pvalues)

# Save regardless of significance
log_mann_whitney_df <- log_mann_whitney_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
dim(log_mann_whitney_df)
print(log_mann_whitney_df[1:5,1:5])
fwrite(log_mann_whitney_df,sprintf("%s/%s_nonPAR_X_eg_correlation_results.tsv", out.dir2, cellLabel),sep="\t",quote=F)

# Save only significant
log_mann_whitney_df_significant <- log_mann_whitney_df[(log_mann_whitney_df$localFDR < 0.05),]
dim(log_mann_whitney_df_significant)
print(log_mann_whitney_df_significant[1:5,1:5])
fwrite(log_mann_whitney_df_significant,sprintf("%s/%s_nonPAR_X_eg_significant_correlation_results.tsv", out.dir2, cellLabel),sep="\t",quote=F)

print('JOB IS DONE!')

quit()

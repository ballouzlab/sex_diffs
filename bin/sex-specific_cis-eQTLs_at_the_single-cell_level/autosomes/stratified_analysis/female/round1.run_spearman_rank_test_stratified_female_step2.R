##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping - round 1
# Author: Seyhan Yazar
# Date: 2022-03
# Description: This R script was written to run an array job per chromosome 
# for all cell types using "round1.run_spearman_rank_test.sh" script
##############################################################################

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-10-12", "female-only")

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
chrNumber <- args[2]
print(chrNumber)

# Example
# cellLabel <- "CD4_CTL"
# chrNumber <- "18"

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)

# Directory paths
data.dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/'

# Input filenames
log_expression_filename <- sprintf("./results/2022-10-12_female-only/%s/%s_log_expression_female.tsv", cellLabel, cellLabel)
covariate_filename <- sprintf("./results/2022-10-12_female-only/%s/%s_covariate_female.tsv", cellLabel, cellLabel)

genotype_filename <- sprintf("%sGenotype_Files/genotype_chr%s.tsv", data.dir, chrNumber)
geneLoc_filename <- sprintf("%sGene_Location_Files/geneloc_chr%s.tsv", data.dir, chrNumber)
snpLoc_filename <- sprintf("%sSNP_Location_Files/snpsloc_chr%s.tsv", data.dir, chrNumber)

# Read in files
## Peer Factors
covariate_df <- fread(covariate_filename)
dim(covariate_df)
print(covariate_df[1:5,1:5])

## Count matrix
log_expression_df <- fread(log_expression_filename)
dim(log_expression_df)
print(log_expression_df[1:5,1:5])

## Genotype file
genotype_df <- fread(genotype_filename)
dim(genotype_df)
print(genotype_df[1:5,1:5])

## Gene location file
geneLoc_df <- fread(geneLoc_filename)
dim(geneLoc_df)
print(geneLoc_df[1:5,1:5])

## SNP location file
snpLoc_df <- fread(snpLoc_filename)
dim(snpLoc_df)
head(snpLoc_df)

# IDs to use
gene_ids <- colnames(log_expression_df[,-"sampleid"])
sample_ids <- log_expression_df$sampleid

# Prepare genotype file
dim(genotype_df)
genotype_df <- genotype_df[(genotype_df$sampleid %in% log_expression_df$sampleid),]
dim(genotype_df)
snp_ids <- colnames(genotype_df[-1])
genotype_df <- genotype_df %>% ungroup()

# Find Gene-SNP pairs for chrNumber
geneLoc_df <- geneLoc_df[(geneLoc_df$geneid %in% gene_ids),]
geneLoc_df$left <- geneLoc_df$start - 1000000
geneLoc_df$right <- geneLoc_df$end + 1000000
gene_snp_df <- geneLoc_df %>% 
  group_by(geneid) %>% 
  group_modify(function(x, y) tibble(snpid = snpLoc_df$snpid[between(as.numeric(snpLoc_df$pos), as.numeric(x["left"]), as.numeric(x["right"]))]), keep=FALSE)
# gene_snp_df <- gene_snp_df %>% ungroup()
gene_snp_df <- inner_join(snpLoc_df, gene_snp_df) %>% select(chr, pos, snpid, geneid)
dim(gene_snp_df)
print(gene_snp_df %>% head)

out.dir = paste0(output,'/', cellLabel)

fwrite(gene_snp_df,sprintf("%s/%s_chr%s_gene_SNP_pairs.tsv", out.dir, cellLabel, chrNumber),sep="\t",quote=F)


# Find residuals for log transformed expressions
calculate_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- log_expression_df %>% select("sampleid",all_of(gene))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,covariate_df,by="sampleid")
  colnames(test_df)[2] <- "expression"

  # Generate model
  model <- lm(expression ~ age + pc1 + pc2 + pc3 + pc4 + pf1 + pf2 , data=test_df)
  residuals=resid(model)
  residuals  
}

log_residual_mat <- sapply(gene_ids,calculate_residuals)
rownames(log_residual_mat) <- log_expression_df$sampleid
log_residual_df <- data.frame(log_residual_mat, check.names=FALSE)
log_residual_df$sampleid <- log_expression_df$sampleid
print(log_residual_df[1:5,1:5])
dim(log_residual_df)

dir.create(paste0(out.dir, "/round1"))

fwrite(log_residual_df,sprintf("%s/round1/%s_chr%s_log_residuals_female.tsv", out.dir, cellLabel, chrNumber),sep="\t",quote=F)

genotype_df <- genotype_df %>% filter(sampleid %in% covariate_df$sampleid)

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
spearman_correlation <- function (x,y) {
  gene <- y$geneid
  # print(gene)
  snp <- y$snpid
  # print(snp)

  # Select values to test
  res_val <- log_residual_df[c("sampleid",gene)]
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

log_spearman_df <- gene_snp_df %>% group_by(geneid,snpid) %>% group_modify(spearman_correlation)

# Calculate the qvalues for pvalues
pvalues <- log_spearman_df$p.value
qobj <- qvalue(p = pvalues)

# Save regardless of significance
log_spearman_df <- log_spearman_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
dim(log_spearman_df)
print(log_spearman_df[1:5,1:5])
fwrite(log_spearman_df,sprintf("%s/round1/%s_chr%s_correlation_results_female.tsv", out.dir, cellLabel,chrNumber),sep="\t",quote=F)

# Save only significant
log_spearman_df_significant <- log_spearman_df[(log_spearman_df$localFDR < 0.05),]
dim(log_spearman_df_significant)
print(log_spearman_df_significant[1:5,1:5])
fwrite(log_spearman_df_significant,sprintf("%s/round1/%s_chr%s_significant_correlation_results_female.tsv", out.dir, cellLabel,chrNumber),sep="\t",quote=F)

print('JOB IS DONE!')

quit()

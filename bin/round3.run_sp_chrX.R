##############################################################################
# Script information                                                      
# Title: Conditional cis-eQTL mapping - round 2
# Author: Seyhan Yazar
# Date: 2021-09-06
# Description: This R script was written to run an array job per chromosome 
# for 26 cell types using "round2.run_spearman_rank_test.sh" script
##############################################################################

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-03-08", "conditional-analysis")

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)
library(future)

# Directory paths
data.dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'

# Input filenames
r1_residual_filename <- sprintf("%s/%s/%s_chrX_eSNP1_adjusted_residuals.tsv",output, cellLabel, cellLabel)
r1_significant_SNPs <- sprintf("%s/%s/%s_chrX_round2_significant_correlation_results.tsv",output,cellLabel, cellLabel)
genotype_filename <- sprintf("data/genotype_chr23.tsv")
geneLoc_filename <- sprintf("data/geneloc_chrX.tsv")
snpLoc_filename <- sprintf("data/onek1k_chr23.bim")

# Read in files
## Count matrix
r1_residual_df <- fread(r1_residual_filename)
print("r1_residaul_df info ...")
dim(r1_residual_df)
print(r1_residual_df[1:5,1:5])

## Significant SNPs file
r1_significant_SNPs_df <- fread(r1_significant_SNPs)
print("r1_significant_SNPs_df info...")
dim(r1_significant_SNPs_df)

## Genotype file
genotype_df <- fread(genotype_filename)
print("genotype_df info...")
dim(genotype_df)

## Gene location file
geneLoc_df <- fread(geneLoc_filename)
dim(geneLoc_df)
colnames(geneLoc_df)[1] <-'geneid'
colnames(geneLoc_df)[3] <-'chr'
print(geneLoc_df[1:5,1:5])

## SNP location file
snpLoc_df <- fread(snpLoc_filename)
dim(snpLoc_df)
colnames(snpLoc_df) <- c('chr', 'snpid','ignore','pos','A1','A2')
head(snpLoc_df)

# Identify the top eSNP for each eGene 
eSNP1 <- r1_significant_SNPs_df %>%
  group_by(geneid) %>%
  arrange(qvalue) %>%
  filter(row_number()==1)

print("eSNP1 info...")
dim(eSNP1)

out.dir = paste0(output,'/', cellLabel)

fwrite(eSNP1,sprintf("%s/%s_chrX_eSNP2.tsv",out.dir,cellLabel),sep="\t",quote=F)

# Remaning eSNPs to test
eSNPs_to_test <- r1_significant_SNPs_df %>% 
    group_by(geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()!=1)  

## Subset residuals for the genes to be tested
dim(r1_residual_df)
sample_ids <- r1_residual_df$sampleid
gene_ids <- eSNP1$geneid
r1_residual_df <- r1_residual_df %>% select (all_of(gene_ids))
r1_residual_df$sampleid <- sample_ids
dim(r1_residual_df)

# Subset genotype file for the significant SNPs
dim(genotype_df)
genotype_df <- genotype_df[(genotype_df$sampleid %in% sample_ids),]
dim(genotype_df)
genotype_df <- genotype_df %>% select(sampleid, r1_significant_SNPs_df$snpid)
genotype_df <- genotype_df %>% ungroup()
snp_ids <- colnames(genotype_df[-1])

# Find residuals after adjustment of lead SNP
calculate_adjusted_residuals <- function (x) {
  gene <- x
  gene <- "HDHD1"

  # select gene to regress
  exprs_val <- r1_residual_df %>% select("sampleid", all_of(gene))

  # select SNP to add
  snp = as.character(eSNP1$snpid[(eSNP1$geneid==gene)])
  snp_genotype = genotype_df %>% select(sampleid, all_of(snp))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,snp_genotype,by="sampleid")
  colnames(test_df)[2] <- "expression"
  colnames(test_df)[3] <- "genotype"
  rownames(test_df) <- test_df$sampleid

  # Generate model
  model <- lm(expression ~ genotype , data=test_df)
  residuals=resid(model)
  residuals
}
options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_residual_mat <- future(sapply(gene_ids,calculate_adjusted_residuals))
adjusted_residual_mat <- value(adjusted_residual_mat)
adjusted_residual_df <- data.frame(adjusted_residual_mat, check.names=FALSE)
adjusted_residual_df$sampleid <- rownames(adjusted_residual_df)
dim(adjusted_residual_df)
fwrite(adjusted_residual_df,sprintf("%s/%s_chrX_eSNP2_adjusted_residuals.tsv",out.dir, cellLabel),sep="\t",quote=F)

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
spearman_correlation <- function (x,y) {
  gene <- y$geneid
  snp <- y$snpid

  # Select values to test
  res_val <- adjusted_residual_df %>% select("sampleid",all_of(gene))
  genotype_val <- genotype_df %>% select("sampleid", all_of(snp))
  
  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by="sampleid")
  colnames(test_df) <- c("sampleid","residual", "SNP")
  
  # Generate model
  model <- cor.test(x=test_df$SNP, y=test_df$residual, method = 'spearman', exact=T)
  model_table <- tidy(model)
  model_table
}

gene_snp_test_df <- eSNPs_to_test %>% select(snpid,geneid)
options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_spearman_df <- future(gene_snp_test_df %>% group_by(snpid,geneid) %>% group_modify(spearman_correlation))
adjusted_spearman_df <- value(adjusted_spearman_df)

# Calculate the qvalues for pvalues
pvalues <- adjusted_spearman_df$p.value
qobj <- qvalue(p = pvalues, pi0=1)

# Save regardless of significance
adjusted_spearman_df <- adjusted_spearman_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
nrow(adjusted_spearman_df)
fwrite(adjusted_spearman_df,sprintf("%s/%s_chrX_round3_correlation_results.tsv",out.dir, cellLabel),sep="\t",quote=F)

# Save only significant
adjusted_spearman_df_significant <- adjusted_spearman_df[(adjusted_spearman_df$localFDR < 0.05),]
nrow(adjusted_spearman_df_significant)
fwrite(adjusted_spearman_df_significant,sprintf("%s/%s_chrX_round3_significant_correlation_results.tsv",out.dir, cellLabel),sep="\t",quote=F)

print('JOB IS DONE!')

quit()
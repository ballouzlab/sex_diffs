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
output <- set_output("2022-09-27", "male-only")

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
chrNumber <- args[2]
print(chrNumber)

# Example
# cellLabel <- "B_intermediate"
# chrNumber <- "22"

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)

# Directory paths
data.dir <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/'
data.dir2 <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'

# Input filenames
expression_filename <- sprintf("%soutputs/filtered_matrices/%s_expression_genes_removed.tsv", data.dir2, cellLabel)
genotype_filename <- sprintf("%sGenotype_Files/genotype_chr%s.tsv", data.dir, chrNumber)
geneLoc_filename <- sprintf("%sGene_Location_Files/geneloc_chr%s.tsv", data.dir, chrNumber)
snpLoc_filename <- sprintf("%sSNP_Location_Files/snpsloc_chr%s.tsv", data.dir, chrNumber)
covariate_filename <- sprintf("%soutputs/peer_factors/%s_peer_factors_after_removing_genes.tsv", data.dir2, cellLabel)

# Read in files
## Peer Factors
covariate_df <- fread(covariate_filename)
dim(covariate_df)
# covariate_df <- covariate_df[,-1]
print(covariate_df[1:5,1:5])

## Count matrix
expression_df <- fread(expression_filename)
dim(expression_df)
expression_df <- expression_df[,-c("V1")]
print(expression_df[1:5,1:5])

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

# Prepare Y variable
## Remove genes with 0 expression in all samples
gene_ids <- colnames(expression_df[,-1])
sample_ids <- expression_df$sampleid

# Work out which ones have colSums != 0:
i <- (colSums(expression_df[,-1], na.rm=T) != 0)
nonzerogenes <- names(i[i==TRUE])
expression_df_nonzero <- expression_df[,nonzerogenes, with = FALSE] %>% 
  mutate (sampleid = sample_ids) %>%
  select(sampleid,everything())
print(expression_df_nonzero[1:5,1:5])
dim(expression_df_nonzero)

# Number of individuals with non-zero expression
numberOfIndividuals <- nrow(expression_df_nonzero)
print(sprintf("Total number of individuals = %s",numberOfIndividuals))
numberOfIndividualsWithNonZeroExpression <- colSums(expression_df_nonzero[,-1]!=0)

# Percentage of individuals with non-zero expression
percentOfIndividualsWithNonZeroExpression <- (numberOfIndividualsWithNonZeroExpression/nrow(expression_df_nonzero))*100

# Mean Expression
meanExpression <-  colMeans(expression_df_nonzero[,-1])

# Expression Variance
expressionVariance <-  apply(expression_df_nonzero[,-1], MARGIN=2, FUN=var, na.rm=TRUE)

# Combine percentage, mean and variance information
extraInfo <- data.frame(colnames(expression_df_nonzero[,-1]), numberOfIndividuals,
  numberOfIndividualsWithNonZeroExpression,
  percentOfIndividualsWithNonZeroExpression,
  meanExpression, expressionVariance)
# fwrite(extraInfo,file=sprintf("%s/round1/%s_%s_extraInfo.tsv", cellLabel, chrNumber),quote=F, sep="\t")

print(sprintf("Number of genes: %s", nrow(extraInfo)))

# Plot gene expression against percent of non-zero individuals
lessthan20 <- extraInfo[(extraInfo$percentOfIndividualsWithNonZeroExpression<20),]
# association <- ggplot(lessthan20, aes(x=percentOfIndividualsWithNonZeroExpression,y=meanExpression)) +
#  geom_point() +
#  theme_classic()
# file=sprintf("%s_%s_expression_lessthan20pctzero.png", cellLabel, chrNumber)
# ggsave(association, here(output,file))

# Filter genes with less than 10 percent individuals with non-zero expression
atleast10percent <- extraInfo[(extraInfo$percentOfIndividualsWithNonZeroExpression>10),]
expression_df_nonzero <- expression_df_nonzero %>% select(rownames(atleast10percent))
expression_df_nonzero$sampleid <- sample_ids
gene_ids <- colnames(expression_df_nonzero)[colnames(expression_df_nonzero)!="sampleid"]

print(sprintf("Number of genes after filtering: %s", nrow(atleast10percent)))

# Prepare genotype file
dim(genotype_df)
genotype_df <- genotype_df[(genotype_df$sampleid %in% expression_df_nonzero$sampleid),]
dim(genotype_df)
snp_ids <- colnames(genotype_df[-1])
genotype_df <- genotype_df %>% ungroup()

# Prepare covariate file
covariate_ids <- colnames(covariate_df[-1])

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
dir.create(out.dir)

fwrite(gene_snp_df,sprintf("%s/%s_chr%s_gene_SNP_pairs.tsv", out.dir, cellLabel, chrNumber),sep="\t",quote=F)

# log+1 transformation
logplusone <- function(x) {log(x[1] + 1)}
log_exprs_mat <- apply(expression_df_nonzero[,1:(length(expression_df_nonzero)-1)],1:2,logplusone)
log_exprs_df <- as.data.frame(log_exprs_mat)
log_exprs_df$sampleid <- sample_ids
print(log_exprs_df[1:5,1:5])

no_males = covariate_df %>% filter(sex==1) %>% nrow() 
no_females = covariate_df %>% filter(sex==2) %>% nrow()
    
covariate_df <- covariate_df %>% filter(sex==1) 
# covariate_df <- covariate_df %>% filter(sex==2)

# if (no_females > no_males) {
#  lm_df_f = sample_n (lm_df_f, no_males)
#   } else {
#      lm_df_f <- lm_df_f
#  }

log_exprs_df <- log_exprs_df %>% filter(sampleid %in% covariate_df$sampleid)

# Find residuals for log transformed expressions
calculate_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- log_exprs_df[,c("sampleid",gene)]

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,covariate_df,by="sampleid")
  colnames(test_df)[2] <- "expression"

  # Generate model
  model <- lm(expression ~ age + pc1 + pc2 + pc3 + pc4 + pf1 + pf2 , data=test_df)
  residuals=resid(model)
  residuals  
}

gene_to_test <- gene_snp_df$geneid %>% unique

log_residual_mat <- sapply(gene_ids,calculate_residuals)
rownames(log_residual_mat) <- covariate_df$sampleid
log_residual_df <- data.frame(log_residual_mat, check.names=FALSE)
log_residual_df$sampleid <- covariate_df$sampleid
print(log_residual_df[1:5,1:5])
dim(log_residual_df)

dir.create(paste0(out.dir, "/round1"))

fwrite(log_residual_df,sprintf("%s/round1/%s_chr%s_log_residuals_male.tsv", out.dir, cellLabel, chrNumber),sep="\t",quote=F)

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
fwrite(log_spearman_df,sprintf("%s/round1/%s_chr%s_correlation_results_male.tsv", out.dir, cellLabel,chrNumber),sep="\t",quote=F)

# Save only significant
log_spearman_df_significant <- log_spearman_df[(log_spearman_df$localFDR < 0.05),]
dim(log_spearman_df_significant)
print(log_spearman_df_significant[1:5,1:5])
fwrite(log_spearman_df_significant,sprintf("%s/round1/%s_chr%s_significant_correlation_results_male.tsv", out.dir, cellLabel,chrNumber),sep="\t",quote=F)

print('JOB IS DONE!')

quit()

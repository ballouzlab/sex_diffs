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

# Example
# cellLabel <- "CD4_Naive"

library(data.table)
library(dplyr)

# celltypes <- c("B_intermediate", "B_memory", "B_naive", 
#    "CD14_Mono", "CD16_Mono", 
#    "CD4_CTL", "CD4_Naive", "CD4_TCM","CD4_TEM",
#    "CD8_Naive", "CD8_TCM", "CD8_TEM", "DC",
#    "dnT","gdT", "MAIT", "Treg",
#    "NK", "NK_CD56bright", "NK_Proliferating", 
#    "Plasmablast"
#    )

# Check the largest cell type
# for (i in seq(1:length(celltypes))) {
#     cellLabel <- celltypes[i]
#     data.dir2 <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'
#     expression_filename <- paste0(data.dir2,"outputs/filtered_matrices/",cellLabel,"_expression_genes_removed.tsv")
#     expression_df <- fread(expression_filename)
#     print(dim(expression_df))  
# }

# Check the ratio of sexes in each cell type
# for (i in seq(1:length(celltypes))) {
#     cellLabel <- celltypes[i]
#     data.dir2 <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'
#     covariate_filename <- sprintf("%soutputs/peer_factors/%s_peer_factors_after_removing_genes.tsv", data.dir2, cellLabel)
#     covariate_df <- fread(covariate_filename)
#     print(table(covariate_df$sex)) 
# }
 
out.dir = paste0(output,'/', cellLabel)
dir.create(out.dir)

data.dir2 <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_June21/'
expression_filename <- sprintf("%soutputs/filtered_matrices/%s_expression_genes_removed.tsv", data.dir2, cellLabel)
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
sample_ids <- expression_df$sampleid

# Individuals with non-zero expression has been eliminated previously
# log+1 transformation
logplusone <- function(x) {log(x[1] + 1)}
log_exprs_mat <- apply(expression_df[,2:(length(expression_df))],1:2,logplusone)
log_exprs_df <- setDT(as.data.frame(log_exprs_mat))
log_exprs_df$sampleid <- sample_ids
print(log_exprs_df[1:5,1:5])

no_males = covariate_df %>% filter(sex==1) %>% nrow() 
no_females = covariate_df %>% filter(sex==2) %>% nrow()
    
covariate_df <- covariate_df %>% filter(sex==2) 
dim(covariate_df)

set.seed(1234)
covariate_df = sample_n (covariate_df, no_males)
dim(covariate_df)

log_exprs_df <- log_exprs_df %>% filter(sampleid %in% covariate_df$sampleid)
fwrite(log_exprs_df,sprintf("%s/%s_log_expression_female.tsv", out.dir, cellLabel),sep="\t",quote=F)

# Make sure the order of log_exprs_df and covariate_df are same
sampleid_df <- log_exprs_df[,"sampleid"]
covariate_df <- left_join(sampleid_df,covariate_df, by="sampleid")
fwrite(covariate_df,sprintf("%s/%s_covariate_female.tsv", out.dir, cellLabel),sep="\t",quote=F)
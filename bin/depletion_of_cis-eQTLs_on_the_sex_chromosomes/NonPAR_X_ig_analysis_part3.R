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

library(broom)
library(ggplot2)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)
library(qvalue)

input.dir = paste0(output,'/NonPAR-X-ig/', cellLabel)

cor_files <- list.files(path=input.dir, pattern=paste0(cellLabel,"_nonPAR_X_ig_correlation_results*"), full.names=T)
cor_list <- lapply(cor_files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
cor_df <- bind_rows(cor_list)

# Calculate the qvalues for pvalues
pvalues <- cor_df$p.value
qobj <- qvalue(p = pvalues)

# Save regardless of significance
log_spearman_df <- cor_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
dim(log_spearman_df)
print(log_spearman_df[1:5,1:5])

# Save only significant
log_spearman_df_significant <- log_spearman_df[(log_spearman_df$localFDR < 0.05),]
dim(log_spearman_df_significant)
print(log_spearman_df_significant[1:5,1:5])
fwrite(log_spearman_df_significant,sprintf("%s/%s_nonPAR_X_ig_significant_correlation_results.tsv", input.dir, cellLabel),sep="\t",quote=F)

print('JOB IS DONE!')

quit()

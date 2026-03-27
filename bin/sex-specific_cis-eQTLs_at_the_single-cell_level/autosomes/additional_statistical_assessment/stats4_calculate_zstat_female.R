#!/usr/bin/env Rscript

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-10-19", "z-testing")
out.dir = paste0(output,"/female")
dir.create(out.dir)

# Call the cell type
args = commandArgs(trailingOnly=TRUE)
celltype1 <- args[1]

# celltype1 <- "DC"

# Upload libraries
library(data.table)
library(dplyr)
library(stringr)

### Analysis for females ###
# Get all the result files (significant eQTls)
female_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only/female_eqtls/female_top_eqtls.tsv"
df_female <- fread(female_file) %>% 
     filter(celltype==celltype1) %>%  
     select(geneid, snpid)

# Retrieve beta and se from matrix eQTL results for both sex

extract_matrix_outputs <- function (x,y) {
     gene_of_interest=y$geneid
     snp_of_interest=y$snpid
     # gene_of_interest="CENPK"
     # snp_of_interest="5:64902091_T"
     chr=sub(":.*","",snp_of_interest)
     
     # female results
     filename_f <- sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-04_matrix-eqtl-analysis-female/%s/%s_chr%s_cis_eqtls.tsv", celltype1, celltype1, chr)
     df_matrix_eqtl_f <- fread(filename_f)
     df_matrix_eqtl_f <-  df_matrix_eqtl_f %>% 
          filter(SNP==snp_of_interest & gene==gene_of_interest) %>% 
          rename(f_beta = beta) %>%
          mutate(f_se=f_beta/abs(`t-stat`)) %>% 
          select(SNP,gene, f_beta, f_se)

     # male results
     filename_m <- sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-04_matrix-eqtl-analysis-male/%s/%s_chr%s_male_cis_eqtls.tsv", celltype1, celltype1, chr)
     df_matrix_eqtl_m <- fread(filename_m)
     df_matrix_eqtl_m <-  df_matrix_eqtl_m %>% 
          filter(SNP==snp_of_interest & gene==gene_of_interest) %>% 
          rename(m_beta = beta) %>%
          mutate(m_se=m_beta/abs(`t-stat`)) %>% 
          select(SNP,gene, m_beta, m_se)
    
     df_matrix_eqtl <- left_join(df_matrix_eqtl_f, df_matrix_eqtl_m, by=c("SNP","gene"))
     df_matrix_eqtl
}

matrix_female_outputs <- df_female %>% group_by(snpid,geneid) %>% group_modify(extract_matrix_outputs)

# Identify the sample size used for celltype1
wd <- '/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-12_female-only'
expression <- fread(sprintf("%s/%s/round1/%s_chr22_log_residuals_female.tsv", wd, celltype1, celltype1))
sample_n <- nrow(expression)

# Add the sample size to the dataframe
matrix_female_outputs$n <- sample_n

# Calculate z-stat for each SNP-gene pair
matrix_female_outputs <-  matrix_female_outputs %>% mutate(z_stat = (f_beta - m_beta) / 
  sqrt(((f_se^2) / n) + ((m_se^2) / n)))

# Calculate the p-values
library(distributions3)
Z <- Normal(0, 1)  # make a standard normal r.v.
matrix_female_outputs <- matrix_female_outputs %>%
     mutate(p_value= 1 - cdf(Z, abs(z_stat)) + cdf(Z, -abs(z_stat)))
matrix_female_outputs <- matrix_female_outputs %>% select(-c("SNP","gene"))

fwrite(matrix_female_outputs, sprintf("%s/%s_female_eQTLs_ztesting_results.tsv", out.dir, celltype1))
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
output <- set_output("2022-12-08", "NK-chrX-by-sex")

library("data.table")
library("dplyr")
library(qvalue)

### Female
f_par1 <- fread("./results/2022-12-09_sex-chr-analysis/PAR1-XX/NK/NK_female_par1_correlation_results.tsv")
f_par1$region <- "XX_PAR1"
f_par2 <- fread("./results/2022-12-09_sex-chr-analysis/PAR2-XX/NK/NK_female_par2_correlation_results.tsv")
f_par2$region <- "XX_PAR2"

nonpar_file_list <- list.files(path="./results/2022-12-09_sex-chr-analysis/NonPAR-XaXi/NK", pattern="NK_nonPAR_XaXi_correlation_results*", full.name=TRUE)
df_nonpar_list <- lapply(nonpar_file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
df_nonpar <- bind_rows(df_nonpar_list)
pvalues <- df_nonpar$p.value
qobj <- qvalue(p = pvalues)
f_nonpar  <- df_nonpar %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
f_nonpar$region <- "XX_nonPAR"

f_nonpar_escp <- fread("./results/2022-12-09_sex-chr-analysis/NonPAR-XaXe/NK/NK_nonPAR_XaXe_correlation_results.tsv")
f_nonpar_escp$region <- "XX_nonPAR_escp"

f_chrx <- bind_rows(f_par1,f_par2, f_nonpar,f_nonpar_escp)

fwrite(f_chrx, file=paste0(output,"/female_chrX_correlations_all_regions_20221209.tsv"), sep="\t")

rm(list=ls())

### Male
m_par1 <- fread("./results/2022-12-09_sex-chr-analysis/PAR1-XY/NK/NK_male_par1_correlation_results.tsv")
m_par1$region <- "XY_PAR1"
m_par2 <- fread("./results/2022-12-09_sex-chr-analysis/PAR2-XY/NK/NK_male_par2_correlation_results.tsv")
m_par2$region <- "XY_PAR2"

nonpar_file_list <- list.files(path="./results/2022-12-09_sex-chr-analysis/NonPAR-X-ig/NK", pattern="NK_nonPAR_X_ig_correlation_results*", full.name=TRUE)
df_nonpar_list <- lapply(nonpar_file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
df_nonpar <- bind_rows(df_nonpar_list)
pvalues <- df_nonpar$p.value
qobj <- qvalue(p = pvalues)
m_nonpar  <- df_nonpar %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
m_nonpar$region <- "XY_nonPAR"

m_chrx <- bind_rows(m_par1,m_par2, m_nonpar)

fwrite(m_chrx, file=paste0(output,"/male_chrX_correlations_all_regions_20221209.tsv"), sep="\t")

f_nonpar <- fread("./results/2022-12-09_sex-chr-analysis/NonPAR-XaXi/NK/NK_nonPAR_XaXi_significant_correlation_results.tsv") %>% 
    group_by(geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1) %>% 
    select(geneid, snpid)

m_nonpar <- fread("./results/2022-12-09_sex-chr-analysis/NonPAR-X-ig/NK/NK_nonPAR_X_ig_significant_correlation_results.tsv") %>%
     group_by(geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1) %>% 
    select(geneid, snpid)

nonpar <- rbind(f_nonpar,m_nonpar) %>% unique()

fwrite(nonpar, "./results/2022-12-08_expr-geno-plots/geneid_snpid_pairs_to_plot.txt", sep=" ")
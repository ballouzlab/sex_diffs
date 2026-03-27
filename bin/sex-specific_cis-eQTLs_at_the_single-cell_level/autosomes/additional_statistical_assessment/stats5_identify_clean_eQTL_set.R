#!/usr/bin/env Rscript

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis")

# Import libraries ------------------------------------
library(tidyverse)
library(dsLib)

# Set output ------------------------------------------
output <- set_output("2022-10-19", "sex-specific-eQTLs-autosomes")

library(data.table)
library(dplyr)

### Analysis for females ###

# z-test results
files <- list.files(path="./results/2022-10-19_z-testing/female", pattern="female_eQTLs_*", full.names=T)
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
celltypes <- gsub("./results/2022-10-19_z-testing/female/*", "", files) %>%
    gsub("*_female_eQTLs_ztesting_results.tsv","",.)
names(dataset) <- celltypes
df_zstat_f <- bind_rows(dataset, .id="celltype")

# pi analysis results
female_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/female_eQTLs_after_pi0_threshold.tsv"
df_pi_f <- fread(female_file) %>%
    select(celltype,snpid, geneid,estimate,p.value,localFDR)

# Combine and filter out using z-test p-value
f_results <- left_join(df_pi_f,df_zstat_f, by=c("celltype", "snpid", "geneid"))
dim(f_results)
f_results <- f_results %>% filter(p_value < 0.05)
dim(f_results)

fwrite(f_results, paste0(output,"/sex_specific_female_eQTLs_not_in_joint_results_and_passed_both_tests.tsv"),sep="\t")

### Analysis for males ###

# z-test results
files <- list.files(path="./results/2022-10-19_z-testing/male", pattern="male_eQTLs_*", full.names=T)
dataset <- lapply(files, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
celltypes <- gsub("./results/2022-10-19_z-testing/male/*", "", files) %>%
    gsub("*_male_eQTLs_ztesting_results.tsv","",.)
names(dataset) <- celltypes
df_zstat_m <- bind_rows(dataset, .id="celltype")

# pi analysis results
male_file <- "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-10-18_pi-analysis/male_eQTLs_after_pi0_threshold.tsv"
df_pi_m <- fread(male_file) %>%
    select(celltype,snpid, geneid,estimate,p.value,localFDR)

# Combine and filter out using z-test p-value
m_results <- left_join(df_pi_m,df_zstat_m, by=c("celltype", "snpid", "geneid"))
dim(m_results)
m_results <- m_results %>% filter(p_value < 0.05)
dim(m_results)

fwrite(m_results, paste0(output,"/sex_specific_male_eQTLs_not_in_joint_results_and_passed_both_tests.tsv"),sep="\t")

#### Additional metrics #####
initial_results <- fread("results/2022-03-14_conditional-analysis-stratified/female_top_eQTLs_26May22.tsv")
initial_results <- initial_results %>% filter(celltype=="CD4_Naive") %>% filter(chr!="23")

initial_results_m <- fread("results/2022-03-14_conditional-analysis-stratified/male_top_eQTLs_26May22.tsv")
initial_results_m <- initial_results_m %>% filter(celltype=="CD4_Naive") %>% filter(chr!="23")

initial_results %>% filter(snpid %in% initial_results_m$snpid) %>% count()
initial_results_m %>% filter(!snpid %in% initial_results$snpid) %>% count()

f_results %>% filter(!snpid %in% shared_eqtls$snpid)

celltype_list <- c("B_intermediate", "B_memory", "B_naive", 
    "CD14_Mono", "CD16_Mono", 
    "CD4_CTL", "CD4_Naive", "CD4_TCM","CD4_TEM",
    "CD8_Naive", "CD8_TCM", "CD8_TEM", "DC",
    #"cDC", "cDC1", "cDC2", "DC", "pDC", "pDC1",
    "dnT","gdT", "MAIT", "Treg",
    "NK", "NK_CD56bright", "NK_Proliferating", 
    "Plasmablast"
    #"Platelet"
    )

count_eqtls <- function (i) {
    # i=celltype_list[1]
    df_f <- f_results %>% filter(celltype==i)
    df_m <- m_results %>% filter(celltype==i)
    f_only = df_f %>% filter(!snpid %in% df_m$snpid) %>% count()
    m_only = df_m %>% filter(!snpid %in% df_f$snpid) %>% count()
    shared = df_f %>% filter(snpid %in% df_m$snpid) %>% count()
    info = data.frame(i,f_only, shared, m_only)
    colnames(info) <- c("celltype","f_only","shared","m_only")
    info
}

bind_rows(lapply(celltype_list, count_eqtls))
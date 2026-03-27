##############################################################################
# Script information                                                      
# Title: Tidy PAR1 region results
# Author: Seyhan Yazar
# Date: 2021-09-01
# Description: 
##############################################################################

library(data.table)
library(dplyr)

# Set working directory
setwd("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/sex_specific_eQTL_analysis/results/2022-12-09_sex-chr-analysis")

# Extract NonPAR XaXe 

file_list <- list.files(path="NonPAR-XaXe", pattern="*_nonPAR_XaXe_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("NonPAR-XaXe/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")

dir.create("NonPAR")

fwrite(df,sprintf("NonPAR/NonPAR-XaXe_all_eqtls.tsv"),sep="\t",quote=F)

eSNP1_XaXe <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

fwrite(eSNP1_XaXe,sprintf("NonPAR/NonPAR-XaXe_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())

# Extract NonPAR XaXi

file_list <- list.files(path="NonPAR-XaXi", pattern="*_nonPAR_XaXi_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("NonPAR-XaXi/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")

fwrite(df,sprintf("NonPAR/NonPAR-XaXi_all_eqtls.tsv"),sep="\t",quote=F)

eSNP1_XaXi <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

fwrite(eSNP1_XaXi,sprintf("NonPAR/NonPAR-XaXi_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())

# Extract NonPAR X_ig

file_list <- list.files(path="NonPAR-X-ig", pattern="*_nonPAR_X_ig_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("NonPAR-X-ig/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")

fwrite(df,sprintf("NonPAR/NonPAR-X-ig_all_eqtls.tsv"),sep="\t",quote=F)

eSNP1_X_ig <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

fwrite(eSNP1_X_ig,sprintf("NonPAR/NonPAR-X-ig_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())

# Extract NonPAR XaXiX

file_list <- list.files(path="NonPAR-XaXiX", pattern="*_nonPAR_XaXiX_significant_correlation_results.tsv", full.name=TRUE, recursive=TRUE)
celltypes <- sub("NonPAR-XaXiX/*", "", file_list) %>% sub("/.*","", .)

df_list <- lapply(file_list, function(x) { tryCatch(fread(x) , error=function(e) NULL)})
names(df_list) <- celltypes
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- bind_rows(df_list, .id="celltype")

fwrite(df,sprintf("NonPAR/NonPAR-XaXiX_all_eqtls.tsv"),sep="\t",quote=F)

eSNP1_X_ig <- df %>%
    group_by(celltype, geneid) %>%
    arrange(qvalue) %>%
    filter(row_number()==1)

fwrite(eSNP1_X_ig,sprintf("NonPAR/NonPAR-XaXiX_top_eqtls.tsv"),sep="\t",quote=F)

rm(list=ls())